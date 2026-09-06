# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""cli.py

Command-line calculator for Mag$\nu$s: computes a single neutrino
oscillation probability (or probability matrix) from the command line,
without writing any Python. Wraps the same ``osc_prob_{2,3,4,5}nu_*``
functions used by the Python API (see :py:mod:`magnus.oscprob`
and :doc:`/cli`), dispatching to the right one based on ``--flavors``,
``--environment``, and ``--scenario``.

Installed as the ``magnus`` console script (see ``[project.scripts]``
in pyproject.toml) and also runnable as ``python -m magnus``.

Routine listings
----------------

    * main - Entry point: parses argv, dispatches, prints the result
    * build_parser - Builds the argparse.ArgumentParser
    * FLAVOR_NAME_TO_INDEX - Maps flavor names (e, mu, tau, s, s1, s2)
           to their globaldefs index
    * ENERGY_UNITS, LENGTH_UNITS - Unit-name to :math:`\text{eV}` / :math:`\text{eV}^{-1}`
           conversion factors
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import argparse
import inspect
import json
import sys
from typing import Optional

import numpy as np

import magnus.oscprob as oscprob
import magnus.globaldefs as gd
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

# Refinement/logging/numerics kwargs that every osc_prob_* wrapper accepts via
# **kwargs even where they are not explicit named parameters (see the "layer
# contract" in docs/source/architecture.rst) -- always safe to forward.
ALWAYS_FORWARD = {'magnus_exp_order', 'n_jobs', 'integration_method', 'rtol', 'atol',
                  'strategy'}


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


def build_parser() -> argparse.ArgumentParser:
    r"""Builds the ``magnus`` command-line argument parser.

    .. versionadded:: 1.0.0

    Returns
    -------
    argparse.ArgumentParser
        The top-level parser, with the ``prob`` subcommand attached.
    """
    parser = argparse.ArgumentParser(
        prog='magnus',
        description="Magνs: neutrino oscillation probabilities via the Magnus expansion.")
    parser.add_argument('-V', '--version', action='version', version=f'magnus {__version__}')
    sub = parser.add_subparsers(dest='command', required=True)

    p = sub.add_parser('prob', help='Compute a single oscillation probability (matrix or channel).')

    g_env = p.add_argument_group('Environment')
    g_env.add_argument('--flavors', type=int, choices=[2, 3, 4, 5], default=3,
        help='Number of neutrino flavors (default: 3).')
    g_env.add_argument('--environment', choices=['vacuum', 'matter', 'earth', 'sun'], default='vacuum',
        help='Propagation environment (default: vacuum).')
    g_env.add_argument('--scenario', choices=['std', 'nsi', 'liv'], default='std',
        help="Physics scenario on top of the environment: 'std' (Standard Model), "
             "'nsi' (non-standard interactions), or 'liv' (Lorentz-invariance violation). "
             "'nsi' is not available with --environment vacuum. Default: std.")
    g_env.add_argument('--density-profile', choices=['constant', 'exp'], default='constant',
        help="Matter density profile, only used with --environment matter: 'constant' "
             "(requires --rho) or 'exp' (requires --rho-central and --l-scale). Default: constant.")
    g_env.add_argument('--nubar', action='store_true',
        help='Compute the probability for antineutrinos instead of neutrinos.')

    g_kin = p.add_argument_group('Energy and baseline')
    g_kin.add_argument('--energy', type=float, required=True, help='Neutrino energy.')
    g_kin.add_argument('--energy-unit', choices=list(ENERGY_UNITS), default='GeV',
        help='Unit of --energy (default: GeV).')
    g_kin.add_argument('--baseline', type=float, default=None,
        help='Baseline / final position. Required for vacuum, matter, and sun, and for earth '
             'when using --costhz. Only computed automatically for earth when both --loc-ini '
             'and --loc-fin are given instead.')
    g_kin.add_argument('--l0', type=float, default=0.0,
        help='Initial position (used by --environment sun and --density-profile exp). '
             'Default: 0.0.')
    g_kin.add_argument('--baseline-unit', choices=list(LENGTH_UNITS), default='km',
        help='Unit of --baseline, --l0, and --l-scale (default: km).')

    g_mat = p.add_argument_group('Matter (--environment matter)')
    g_mat.add_argument('--rho', type=float, default=None,
        help='Matter density (constant profile).')
    g_mat.add_argument('--rho-central', type=float, default=None,
        help='Matter density at the center of the profile, l=0 (exponential profile).')
    g_mat.add_argument('--l-scale', type=float, default=None,
        help='Length scale of the exponential density decrease (exponential profile).')
    g_mat.add_argument('--density-unit', choices=['g/cm3', 'natural'], default='g/cm3',
        help='Unit of --rho/--rho-central: g/cm3 (converted internally) or natural units '
             '(eV^4). Default: g/cm3.')
    g_mat.add_argument('--ratio-n-to-p', type=float, default=1.0,
        help='Ratio of the number of neutrons to protons in matter. Default: 1.0.')
    g_mat.add_argument('--electron-fraction', type=float, default=0.5,
        help='Electron fraction of matter. Default: 0.5.')

    g_earth = p.add_argument_group('Earth (--environment earth)')
    g_earth.add_argument('--costhz', type=float, default=None,
        help='Cosine of the neutrino zenith angle.')
    g_earth.add_argument('--loc-ini', default=None,
        help='Initial location name (e.g. fermilab); see magnus.earth.loc_coords_dms. '
             'Must be given together with --loc-fin, as an alternative to --costhz.')
    g_earth.add_argument('--loc-fin', default=None,
        help='Final location name; see --loc-ini.')

    g_osc = p.add_argument_group('Standard oscillation parameters (2-flavor)')
    g_osc.add_argument('--angles', default='sin', choices=list(gd.ANGLE_CONVENTIONS),
        help='Convention for every mixing angle below (--sth, --s12.., --sxi..): sin (default) '
             'their sines, sin2 their sines squared -- the form global fits report -- rad the '
             'angles in radians, deg in degrees. Under deg the CP phases (--dcp, --d14, ...) '
             'are read as degrees too; otherwise they stay in radians.')
    g_osc.add_argument('--sth', type=float, default=None,
        help='Mixing angle theta, in the convention set by --angles (required for --flavors 2).')
    g_osc.add_argument('--dm2', type=float, default=None, dest='Dm2',
        help='Mass-squared difference Delta m^2 (required for --flavors 2).')

    g_osc3 = p.add_argument_group('Standard oscillation parameters (3+ flavors)')
    g_osc3.add_argument('--s12', type=float, default=None, help='Mixing angle theta_12, per --angles. Default: NuFit 6.1.')
    g_osc3.add_argument('--s23', type=float, default=None, help='Mixing angle theta_23, per --angles. Default: NuFit 6.1.')
    g_osc3.add_argument('--s13', type=float, default=None, help='Mixing angle theta_13, per --angles. Default: NuFit 6.1.')
    g_osc3.add_argument('--dcp', type=float, default=None, dest='dCP',
        help='delta_CP [radian, or degree with --angles deg]. Default: NuFit 6.1.')
    g_osc3.add_argument('--dm21', type=float, default=None, dest='D21',
        help='Mass-squared difference Delta m^2_21. Default: NuFit 6.1.')
    g_osc3.add_argument('--dm31', type=float, default=None, dest='D31',
        help='Mass-squared difference Delta m^2_31. Default: NuFit 6.1.')
    g_osc3.add_argument('--osc-params-set', default='OSC_PARAMS_DEFAULT',
        dest='default_osc_params_set_name',
        choices=sorted(gd.OSC_PARAMS_PREDEFINED),
        help='Predefined set used to fill in any of s12/s23/s13/dCP/D21/D31 left unspecified: '
             'normal ordering (..._NO) or inverted ordering (..._IO).  OSC_PARAMS_DEFAULT is '
             'NuFit 6.1 NO.  Taken from globaldefs.OSC_PARAMS_PREDEFINED rather than listed here, '
             'because a hand-written list went stale: it offered only the 6.0 sets, so asking for '
             'inverted ordering silently dropped a release behind the default.')

    g_osc4 = p.add_argument_group('Additional sterile mixing (4+ flavors)')
    g_osc4.add_argument('--s14', type=float, default=0.0, help='Mixing angle theta_14, per --angles. Default: 0.0.')
    g_osc4.add_argument('--d14', type=float, default=0.0, help='delta_14 [radian, or degree with --angles deg]. Default: 0.0.')
    g_osc4.add_argument('--s24', type=float, default=0.0, help='Mixing angle theta_24, per --angles. Default: 0.0.')
    g_osc4.add_argument('--d24', type=float, default=0.0, help='delta_24 [radian, or degree with --angles deg]. Default: 0.0.')
    g_osc4.add_argument('--s34', type=float, default=0.0, help='Mixing angle theta_34, per --angles. Default: 0.0.')
    g_osc4.add_argument('--dm41', type=float, default=0.0, dest='D41',
        help='Mass-squared difference Delta m^2_41. Default: 0.0.')

    g_osc5 = p.add_argument_group('Additional sterile mixing (5 flavors)')
    g_osc5.add_argument('--s15', type=float, default=0.0, help='Mixing angle theta_15, per --angles. Default: 0.0.')
    g_osc5.add_argument('--d15', type=float, default=0.0, help='delta_15 [radian, or degree with --angles deg]. Default: 0.0.')
    g_osc5.add_argument('--s25', type=float, default=0.0, help='Mixing angle theta_25, per --angles. Default: 0.0.')
    g_osc5.add_argument('--s35', type=float, default=0.0, help='Mixing angle theta_35, per --angles. Default: 0.0.')
    g_osc5.add_argument('--d35', type=float, default=0.0, help='delta_35 [radian, or degree with --angles deg]. Default: 0.0.')
    g_osc5.add_argument('--dm51', type=float, default=0.0, dest='D51',
        help='Mass-squared difference Delta m^2_51. Default: 0.0.')

    g_nsi = p.add_argument_group('NSI parameters (--scenario nsi)')
    g_nsi.add_argument('--eps-aa', type=float, default=0.0, help='2-flavor diagonal NSI coupling.')
    g_nsi.add_argument('--eps-ab', type=float, default=0.0, help='2-flavor off-diagonal NSI coupling.')
    g_nsi.add_argument('--eps-ee', type=float, default=0.0, help='Diagonal NSI coupling of nu_e.')
    g_nsi.add_argument('--eps-em', type=float, default=0.0, help='Off-diagonal (e-mu) NSI coupling.')
    g_nsi.add_argument('--eps-et', type=float, default=0.0, help='Off-diagonal (e-tau) NSI coupling.')
    g_nsi.add_argument('--eps-mm', type=float, default=0.0, help='Diagonal NSI coupling of nu_mu.')
    g_nsi.add_argument('--eps-mt', type=float, default=0.0, help='Off-diagonal (mu-tau) NSI coupling.')
    g_nsi.add_argument('--eps-tt', type=float, default=0.0, help='Diagonal NSI coupling of nu_tau.')
    g_nsi.add_argument('--eps-es', type=float, default=0.0, help='(4nu) Off-diagonal (e-s) NSI coupling.')
    g_nsi.add_argument('--eps-ms', type=float, default=0.0, help='(4nu) Off-diagonal (mu-s) NSI coupling.')
    g_nsi.add_argument('--eps-ts', type=float, default=0.0, help='(4nu) Off-diagonal (tau-s) NSI coupling.')
    g_nsi.add_argument('--eps-ss', type=float, default=0.0, help='(4nu) Diagonal NSI coupling of nu_s.')
    g_nsi.add_argument('--eps-es1', type=float, default=0.0, help='(5nu) Off-diagonal (e-s1) NSI coupling.')
    g_nsi.add_argument('--eps-es2', type=float, default=0.0, help='(5nu) Off-diagonal (e-s2) NSI coupling.')
    g_nsi.add_argument('--eps-ms1', type=float, default=0.0, help='(5nu) Off-diagonal (mu-s1) NSI coupling.')
    g_nsi.add_argument('--eps-ms2', type=float, default=0.0, help='(5nu) Off-diagonal (mu-s2) NSI coupling.')
    g_nsi.add_argument('--eps-ts1', type=float, default=0.0, help='(5nu) Off-diagonal (tau-s1) NSI coupling.')
    g_nsi.add_argument('--eps-ts2', type=float, default=0.0, help='(5nu) Off-diagonal (tau-s2) NSI coupling.')
    g_nsi.add_argument('--eps-s1s1', type=float, default=0.0, help='(5nu) Diagonal NSI coupling of nu_s1.')
    g_nsi.add_argument('--eps-s1s2', type=float, default=0.0, help='(5nu) Off-diagonal (s1-s2) NSI coupling.')
    g_nsi.add_argument('--eps-s2s2', type=float, default=0.0, help='(5nu) Diagonal NSI coupling of nu_s2.')

    g_liv = p.add_argument_group('LIV parameters (--scenario liv)')
    g_liv.add_argument('--sxi', type=float, default=0.0, help='2-flavor LIV mixing angle xi, per --angles.')
    g_liv.add_argument('--sxi12', type=float, default=0.0, help='LIV mixing angle xi_12, per --angles.')
    g_liv.add_argument('--sxi23', type=float, default=0.0, help='LIV mixing angle xi_23, per --angles.')
    g_liv.add_argument('--sxi13', type=float, default=0.0, help='LIV mixing angle xi_13, per --angles.')
    g_liv.add_argument('--dxicp', type=float, default=0.0, dest='dxiCP',
        help='(3nu) LIV CP-violation phase [radian].')
    g_liv.add_argument('--dxi13', type=float, default=0.0,
        help='(4/5nu) LIV CP-violation phase [radian] (replaces --dxicp).')
    g_liv.add_argument('--sxi14', type=float, default=0.0, help='(4/5nu) LIV mixing angle xi_14, per --angles.')
    g_liv.add_argument('--dxi14', type=float, default=0.0, help='(4/5nu) LIV CP-violation phase [radian].')
    g_liv.add_argument('--sxi24', type=float, default=0.0, help='(4/5nu) LIV mixing angle xi_24, per --angles.')
    g_liv.add_argument('--dxi24', type=float, default=0.0, help='(4/5nu) LIV CP-violation phase [radian].')
    g_liv.add_argument('--sxi34', type=float, default=0.0, help='(4/5nu) LIV mixing angle xi_34, per --angles.')
    g_liv.add_argument('--sxi15', type=float, default=0.0, help='(5nu) LIV mixing angle xi_15, per --angles.')
    g_liv.add_argument('--dxi15', type=float, default=0.0, help='(5nu) LIV CP-violation phase [radian].')
    g_liv.add_argument('--sxi25', type=float, default=0.0, help='(5nu) LIV mixing angle xi_25, per --angles.')
    g_liv.add_argument('--sxi35', type=float, default=0.0, help='(5nu) LIV mixing angle xi_35, per --angles.')
    g_liv.add_argument('--dxi35', type=float, default=0.0, help='(5nu) LIV CP-violation phase [radian].')
    g_liv.add_argument('--b1', type=float, default=0.0, help='LIV eigenvalue b1.')
    g_liv.add_argument('--b2', type=float, default=0.0, help='LIV eigenvalue b2.')
    g_liv.add_argument('--b3', type=float, default=0.0, help='LIV eigenvalue b3.')
    g_liv.add_argument('--b4', type=float, default=0.0, help='LIV eigenvalue b4.')
    g_liv.add_argument('--b5', type=float, default=0.0, help='LIV eigenvalue b5.')
    g_liv.add_argument('--liv-lambda', type=float, default=1.0, dest='Lambda',
        help='LIV energy scale Lambda. Default: 1.0.')
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
        help='Highest order of the Magnus expansion (1-8). Default: 4.')
    g_num.add_argument('--integration-method', choices=['gl', 'trapezoid', 'simpson'], default='gl',
        help="Quadrature method. 'gl' (Gauss-Legendre collocation) needs only 1-4 Hamiltonian "
             "evaluations per slab and matches its quadrature order to the expansion order, so "
             "it is both the fastest and the most accurate for a smooth Hamiltonian. "
             "'trapezoid'/'simpson' sample a uniform grid of --n-tpts-per-slab points instead, "
             "and are the safer choice if the Hamiltonian is not smooth within a slab. "
             "Default: gl.")
    g_num.add_argument('--rtol', type=float, default=1.e-3,
        help='Relative tolerance on the agreement between successive refinement levels -- a stopping rule, not a guaranteed accuracy. Default: 1e-3.')
    g_num.add_argument('--atol', type=float, default=1.e-3,
        help='Absolute tolerance on the same agreement; see --rtol. Default: 1e-3.')
    g_num.add_argument('--n-jobs', type=int, default=1, dest='n_jobs',
        help='Number of parallel joblib workers. Default: 1.')
    g_num.add_argument('--strategy', choices=['auto', 'hybrid', 'magnus'], default='auto',
        help="How to propagate a position-dependent Hamiltonian: 'magnus' uses only the "
             "Magnus-expansion machinery; 'hybrid' also tries adiabatic transport with a "
             "Magnus patch at each non-adiabatic window, warning if it cannot certify the "
             "result; 'auto' tries hybrid and falls back to magnus silently. Ignored for "
             "vacuum and constant-density environments. Default: auto.")
    g_num.add_argument('--verbose', type=int, default=0, choices=[0, 1, 2],
        help='Verbosity level. Default: 0.')

    g_out = p.add_argument_group('Output')
    g_out.add_argument('--json', action='store_true', help='Print the result as JSON instead of a table.')
    g_out.add_argument('--precision', type=int, default=4, help='Decimal digits shown in table output. Default: 4.')

    return parser


def _std_osc_kwargs(flavors: int, args: argparse.Namespace) -> dict:
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
    if flavors == 2:
        return {'sxi': args.sxi, 'b1': args.b1, 'b2': args.b2, 'Lambda': args.Lambda,
                'n_liv': args.n_liv}
    kw = {'sxi12': args.sxi12, 'sxi23': args.sxi23, 'sxi13': args.sxi13,
          'b1': args.b1, 'b2': args.b2, 'b3': args.b3, 'Lambda': args.Lambda, 'n_liv': args.n_liv}
    if flavors == 3:
        kw['dxiCP'] = args.dxiCP
    if flavors >= 4:
        kw.update({'dxi13': args.dxi13, 'sxi14': args.sxi14, 'dxi14': args.dxi14,
                   'sxi24': args.sxi24, 'dxi24': args.dxi24, 'sxi34': args.sxi34, 'b4': args.b4})
    if flavors == 5:
        kw.update({'sxi15': args.sxi15, 'dxi15': args.dxi15, 'sxi25': args.sxi25,
                   'sxi35': args.sxi35, 'dxi35': args.dxi35, 'b5': args.b5})
    return kw


def _env_kwargs(environment: str, density_profile: str, args: argparse.Namespace,
                 baseline_ev: Optional[float], l0_ev: float) -> dict:
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
        if not using_locations:
            if args.costhz is None:
                raise SystemExit("magnus prob: --environment earth requires either --costhz "
                                  "(together with --baseline) or both --loc-ini and --loc-fin.")
            if baseline_ev is None:
                raise SystemExit("magnus prob: --baseline is required together with --costhz "
                                  "(only --loc-ini/--loc-fin compute the baseline automatically).")
        return {'costhz': args.costhz, 'loc_ini': args.loc_ini, 'loc_fin': args.loc_fin,
                'L': baseline_ev}
    if environment == 'sun':
        if baseline_ev is None:  # pragma: no cover - pre-empted, see below
            # Unreachable from the command line: main() rejects a missing --baseline for
            # every environment except earth before it calls this function, so a solar run
            # without one has already exited.  Kept because this function is the one that
            # knows what the sun branch needs, and a future caller reaching it by another
            # route should still get a clear message rather than a KeyError downstream.
            raise SystemExit("magnus prob: --baseline is required for --environment sun.")
        return {'L': baseline_ev, 'L0': l0_ev}
    raise AssertionError(environment)  # pragma: no cover


def _wrapper_name(flavors: int, environment: str, scenario: str, density_profile: str) -> str:
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
    universally-forwarded refinement/logging/numerics kwargs)."""
    sig = inspect.signature(fn)
    explicit_names = {n for n, par in sig.parameters.items()
                       if par.kind != inspect.Parameter.VAR_KEYWORD}
    kwargs = {k: v for k, v in candidate_kwargs.items()
              if (k in ALWAYS_FORWARD) or (k in explicit_names)}
    return fn(**kwargs)


def _format_table(P: np.ndarray, flavors: int, precision: int) -> str:
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


def main(argv=None) -> int:
    r"""Entry point for the ``magnus`` console script / ``python -m magnus``.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    argv : list of str, optional
        Arguments to parse instead of ``sys.argv[1:]`` (mainly for testing).

    Returns
    -------
    int
        Process exit code (0 on success).
    """
    parser = build_parser()
    args = parser.parse_args(argv)

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

    fn_name = _wrapper_name(flavors, environment, scenario, args.density_profile)
    fn = getattr(oscprob, fn_name)

    candidate = {'energy': energy_ev}
    candidate.update(_env_kwargs(environment, args.density_profile, args, baseline_ev, l0_ev))
    candidate.update(_std_osc_kwargs(flavors, args))
    if scenario == 'nsi':
        candidate.update(_nsi_kwargs(flavors, args))
    elif scenario == 'liv':
        candidate.update(_liv_kwargs(flavors, args))
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
                                           and args.density_profile == 'exp'):
        candidate['strategy'] = args.strategy

    try:
        P = _call(fn, candidate)
    except ValueError as error:
        # The library validates its own inputs and raises; surface that as a clean CLI error
        # (exit code 2, like any other argument problem) rather than a raw traceback.
        parser.error(str(error))

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
