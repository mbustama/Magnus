# -*- coding: utf-8 -*-
r"""__init__.py

Contains physical constants and unit-conversion constants.

This module contains values of physical constants and unit-conversion
factors used by the various modules of Magnus: unit conversions (km,
cm, GeV, etc., to natural units of eV), fundamental constants (G_F,
particle masses, Avogadro's number), Earth/Sun radii and reference
densities, flavor index constants (NUE, NUMU, NUTAU, NUS), predefined
oscillation/NSI/LIV parameter sets (e.g., NuFit 6.0), and ANSI terminal
color codes (class ``cstyle``) used to format warning/error messages.

Routine listings
----------------

    * cstyle - ANSI terminal color-code constants

The remaining module-level names are physical constants and
unit-conversion factors, not routines; see the module source for the
full list.
"""


__version__ = "1.0"
__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


# from numpy import *
import numpy as np

import os
import platform

# If on Windows, need to call os.system() to print in color in stdout
if platform.system() == 'Windows':
    os.system("")

# Class of different styles
class cstyle():
    r"""ANSI escape-code constants for colored/styled terminal output.

    Used to format the warning/error/tolerance messages printed by
    ``oscprob.py`` (e.g., ``gd.WARNING_MSG_IN_COLOR``,
    ``gd.ERROR_MSG_IN_COLOR``). Has no effect on Windows terminals unless
    ``os.system("")`` has been called first, which this module does at
    import time.

    .. versionadded:: 0.10.0
    """

    CEND      = '\33[0m'
    CBOLD     = '\33[1m'
    CITALIC   = '\33[3m'
    CURL      = '\33[4m'
    CBLINK    = '\33[5m'
    CBLINK2   = '\33[6m'
    CSELECTED = '\33[7m'

    CBLACK  = '\33[30m'
    CRED    = '\33[31m'
    CGREEN  = '\33[32m'
    CYELLOW = '\33[33m'
    CBLUE   = '\33[34m'
    CVIOLET = '\33[35m'
    CBEIGE  = '\33[36m'
    CWHITE  = '\33[37m'

    CBLACKBG  = '\33[40m'
    CREDBG    = '\33[41m'
    CGREENBG  = '\33[42m'
    CYELLOWBG = '\33[43m'
    CBLUEBG   = '\33[44m'
    CVIOLETBG = '\33[45m'
    CBEIGEBG  = '\33[46m'
    CWHITEBG  = '\33[47m'

    CGREY    = '\33[90m'
    CRED2    = '\33[91m'
    CGREEN2  = '\33[92m'
    CYELLOW2 = '\33[93m'
    CBLUE2   = '\33[94m'
    CVIOLET2 = '\33[95m'
    CBEIGE2  = '\33[96m'
    CWHITE2  = '\33[97m'

    CGREYBG    = '\33[100m'
    CREDBG2    = '\33[101m'
    CGREENBG2  = '\33[102m'
    CYELLOWBG2 = '\33[103m'
    CBLUEBG2   = '\33[104m'
    CVIOLETBG2 = '\33[105m'
    CBEIGEBG2  = '\33[106m'
    CWHITEBG2  = '\33[107m'
    # BLACK = '\033[30m'
    # RED = '\033[31m'
    # GREEN = '\033[32m'
    # YELLOW = '\033[33m'
    # BLUE = '\033[34m'
    # MAGENTA = '\033[35m'
    # CYAN = '\033[36m'
    # WHITE = '\033[37m'
    # UNDERLINE = '\033[4m'
    # RESET = '\033[0m'


WARNING_MSG_NO_COLOR = "Warning:"

WARNING_MSG_IN_COLOR = cstyle.CVIOLETBG + "Warning:" + cstyle.CEND

ERROR_MSG_NO_COLOR = "Error in magnus:"

ERROR_MSG_IN_COLOR = cstyle.CREDBG + "Error in magnus:" + cstyle.CEND

TOL_MSG_NO_COLOR = "Requested tolerance achieved"

TOL_MSG_IN_COLOR = cstyle.CGREENBG + "Requested tolerance achieved" + cstyle.CEND


MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = 5
r"""float: Module-level constant

Maximum number of flavors for which we have hard-coded routines in the oscprob module.
Units: [Adimensional]
"""


MAGNUS_EXP_ORDER_MAX = 6
r"""float: Module-level constant

Maximum order of the Magnus expansion currently supported.
Units: [Adimensional]
"""


CONV_KM_TO_INV_EV = 5.06773e9
r"""float: Module-level constant

Multiplicative conversion factor from km to eV^{-1}.
Units: [km^{-1} eV^{-1}].
"""

UNIT_KM = CONV_KM_TO_INV_EV
r"""float: Module-level constant

Multiplicative conversion factor from km to eV^{-1}.  Alias for CONV_KM_TO_INV_EV.
Units: [km^{-1} eV^{-1}].
"""

CONV_CM_TO_INV_EV = CONV_KM_TO_INV_EV*1.e-5
r"""float: Module-level constant

Multiplicative conversion factor from cm to eV^{-1}.
Units: [cm^{-1} eV^{-1}]
"""

UNIT_CM = CONV_CM_TO_INV_EV
r"""float: Module-level constant

Multiplicative conversion factor from cm to eV^{-1}.  Alias for CONV_CM_TO_INV_EV.
Units: [cm^{-1} eV^{-1}]
"""

CONV_CM3_TO_INV_EV3 = np.power(CONV_CM_TO_INV_EV, 3.0)
r"""float: Module-level constant

Multiplicative conversion factor from cm^3 to eV^{-3}.
Units: [cm^{-3} eV^{-3}]
"""

UNIT_CM3 = CONV_CM3_TO_INV_EV3
r"""float: Module-level constant

Multiplicative conversion factor from cm^3 to eV^{-3}.  Alias for CONV_CM3_TO_INV_EV3.
Units: [cm^{-3} eV^{-3}]
"""

CONV_INV_EV_TO_CM = 1./CONV_CM_TO_INV_EV
r"""float: Module-level constant

Multiplicative conversion factor from eV^{-1} to cm.
Units: [eV cm]
"""

UNIT_PER_CM3 = 1.0/CONV_CM3_TO_INV_EV3
r"""float: Module-level constant

Multiplicative conversion factor from cm^{-3} to eV^3. 
Units: [cm^3 eV^3]
"""

CONV_EV_TO_G = 1.783e-33
r"""float: Module-level constant

Multiplicative conversion factor from eV^{-1} to grams.
Units: [g eV^{-1}]
"""

CONV_G_TO_EV = 1./CONV_EV_TO_G
r"""float: Module-level constant

Multiplicative conversion factor from grams to eV^{-1}.
Units: [eV g^{-1}]
"""

UNIT_G_PER_CM3 = CONV_G_TO_EV/CONV_CM3_TO_INV_EV3
r"""float: Module-level constant

Multiplicative conversion factor from g cm^{-3} to eV^4.
Units: [g^{-1} cm^3 eV^4]
"""


SQRT_OF_2 = np.sqrt(2.0)
r"""float: Module-level constant

Square root of 2..
Units: [Adimensional]
"""

GF = 1.1663787e-23
r"""float: Module-level constant

Fermi constant.
Units: [eV^{-2}]
"""

MASS_ELECTRON = 0.5109989461e6
r"""float: Module-level constant

Electron mass.
Units: [eV]
"""

MASS_PROTON = 938.272046e6
r"""float: Module-level constant

Proton mass.
Units: [eV]
"""

MASS_NEUTRON = 939.565379e6
r"""float: Module-level constant

Neutron mass.
Units: [eV]
"""

ELECTRON_FRACTION_EARTH_CRUST = 0.5
r"""float: Module-level constant

Electron fraction in the Earth's crust.
Units: [Adimensional]
"""

DENSITY_MATTER_CRUST_G_PER_CM3 = 3.0
r"""float: Module-level constant

Average matter density in the Earth's crust.
Units: [g cm^{-3}]
"""

N_AV = 6.02214076e23
r"""float: Module-level constant

Avogadro constant
Units: [mol^{-1}]
"""

# NUM_DENSITY_E_EARTH_CRUST = DENSITY_MATTER_CRUST_G_PER_CM3 * CONV_G_TO_EV \
#                             / ((MASS_PROTON+MASS_NEUTRON)/2.0) \
#                             * ELECTRON_FRACTION_EARTH_CRUST \
#                             / pow(CONV_CM_TO_INV_EV, 3.0)
NUM_DENSITY_E_EARTH_CRUST = DENSITY_MATTER_CRUST_G_PER_CM3 * CONV_G_TO_EV \
                            / ((MASS_PROTON+MASS_NEUTRON)/2.0) \
                            * ELECTRON_FRACTION_EARTH_CRUST \
                            / pow(CONV_CM_TO_INV_EV, 3.0)
r"""float: Module-level constant

Electron number density in the Earth's crust
Units: [eV^3]
"""

VCC_EARTH_CRUST = np.sqrt(2.0)*GF*NUM_DENSITY_E_EARTH_CRUST
r"""float: Module-level constant

Charged-current matter potential in the Earth's crust.
Units: [eV]
"""

EARTH_RADIUS = 6371.0
r"""float: Module-level constant

Average Earth radius.
Units: [km]
"""

SUN_RADIUS = 6.947e5
r"""float: Module-level constant

Average solar radius.
Units: [km]
"""

NUM_DENSITY_E_SUN_CENTRAL = 245.0*N_AV*UNIT_PER_CM3
r"""float: Module-level constant

Electron number density at the center of the Sun.
Units: [eV^3]
"""

L_SCALE_SUN = SUN_RADIUS/10.54*UNIT_KM
r"""float: Module-level constant

Electron number density at the center of the Sun.
Units: [eV^{-1}]
"""

NUE = 0
r"""float: Module-level constant

Index used to denote nu_e flavor when computing probabilities.
Units: [Adimensional]
"""

NUMU = 1
r"""float: Module-level constant

Index used to denote nu_mu flavor when computing probabilities.
Units: [Adimensional]
"""

NUTAU = 2
r"""float: Module-level constant

Index used to denote nu_tau flavor when computing probabilities.
Units: [Adimensional]
"""

NUS = 3
r"""float: Module-level constant

Index used to denote the sterile flavor in when computing four-neutrino
(3+1) probabilities.
Units: [Adimensional]
"""

NUS1 = 3
r"""float: Module-level constant

Index used to denote the first sterile flavor in when computing 
five-neutrino (3+2) probabilities.
Units: [Adimensional]
"""

NUS2 = 4
r"""float: Module-level constant

Index used to denote the second sterile flavor in when computing 
five-neutrino (3+2) probabilities.
Units: [Adimensional]
"""


UNIT_KEV = 1.e3
UNIT_MEV = 1.e6
UNIT_GEV = 1.e9
UNIT_TEV = 1.e12
UNIT_PEV = 1.e15
UNIT_EEV = 1.e18


S12_NO_BF_NUFIT_6_0 = np.sqrt(0.308)
r"""float: Module-level constant

Lepton mixing angle :math:`\sin\theta_{12}`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [Adimensional]
"""

S23_NO_BF_NUFIT_6_0 = np.sqrt(0.470)
r"""float: Module-level constant

Lepton mixing angle :math:`\sin\theta_{23}`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [Adimensional]
"""

S13_NO_BF_NUFIT_6_0 = np.sqrt(2.215e-2)
r"""float: Module-level constant

Lepton mixing angle :math:`\sin\theta_{13}`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [Adimensional]
"""

DCP_NO_BF_NUFIT_6_0 = 212./180.*np.pi
r"""float: Module-level constant

Lepton CP-violation phase :math:`\delta_\text{CP}`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [radian]
"""

D21_NO_BF_NUFIT_6_0 = 7.49e-5
r"""float: Module-level constant

Mass-squared difference :math:`\Delta m_{21}^2`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

D31_NO_BF_NUFIT_6_0 = 2.513e-3
r"""float: Module-level constant

Mass-squared difference :math:`\Delta m_{31}^2`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

S12_IO_BF_NUFIT_6_0 = np.sqrt(0.308)
r"""float: Module-level constant

Lepton mixing angle :math:`\sin\theta_{12}`, best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [Adimensional]
"""

S23_IO_BF_NUFIT_6_0 = np.sqrt(0.550)
r"""float: Module-level constant

Lepton mixing angle :math:`\sin\theta_{23}`, best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [Adimensional]
"""

S13_IO_BF_NUFIT_6_0 = np.sqrt(2.231e-2)
r"""float: Module-level constant

Lepton mixing angle :math:`\sin\theta_{13}`, best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [Adimensional]
"""

DCP_IO_BF_NUFIT_6_0 = 274./180.*np.pi
r"""float: Module-level constant

Lepton CP-violation phase :math:`\delta_\text{CP}`, best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [radian]
"""

D21_IO_BF_NUFIT_6_0 = 7.49e-5
r"""float: Module-level constant

Mass-squared difference :math:`\Delta m_{21}^2`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

D32_IO_BF_NUFIT_6_0 = -2.484e-3
r"""float: Module-level constant

Mass-squared difference :math:`\Delta m_{32}^2`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

D31_IO_BF_NUFIT_6_0 = D32_IO_BF_NUFIT_6_0+D21_IO_BF_NUFIT_6_0
r"""float: Module-level constant

Mass-squared difference :math:`\Delta m_{31}^2`, best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [eV^2]
"""

OSC_PARAMS_NU_FIT_6_0_SK_NO = {
    'name': 'OSC_PARAMS_NU_FIT_6_0_NO',
    'description': 'NuFit 6.0, NO, with SK atmospheric data',
    's12': S12_NO_BF_NUFIT_6_0,
    's23': S23_NO_BF_NUFIT_6_0,
    's13': S13_NO_BF_NUFIT_6_0,
    'dCP': DCP_NO_BF_NUFIT_6_0,
    'D21': D21_NO_BF_NUFIT_6_0,
    'D31': D31_NO_BF_NUFIT_6_0
}

OSC_PARAMS_NU_FIT_6_0_SK_IO = {
    'name': 'OSC_PARAMS_NU_FIT_6_0_IO',
    'description': 'NuFit 6.0, IO, with SK atmospheric data',
    's12': S12_IO_BF_NUFIT_6_0,
    's23': S23_IO_BF_NUFIT_6_0,
    's13': S13_IO_BF_NUFIT_6_0,
    'dCP': DCP_IO_BF_NUFIT_6_0,
    'D21': D21_IO_BF_NUFIT_6_0,
    'D31': D31_IO_BF_NUFIT_6_0
}

OSC_PARAMS_PREDEFINED = {
    'OSC_PARAMS_DEFAULT': OSC_PARAMS_NU_FIT_6_0_SK_NO,
    'OSC_PARAMS_NU_FIT_6_0_SK_NO': OSC_PARAMS_NU_FIT_6_0_SK_NO,
    'OSC_PARAMS_NU_FIT_6_0_SK_IO': OSC_PARAMS_NU_FIT_6_0_SK_IO
}


EPS_EE = 0.06
r"""float: Module-level constant

Total NSI strength parameter eps_ee computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_EM = -0.06
r"""float: Module-level constant

Total NSI strength parameter eps_em computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_ET = 0.0
r"""float: Module-level constant

Total NSI strength parameter eps_et computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_MM = 1.2
r"""float: Module-level constant

Total NSI strength parameter eps_mm computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_MT = 0.0
r"""float: Module-level constant

Total NSI strength parameter eps_mt computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_TT = 0.0
r"""float: Module-level constant

Total NSI strength parameter eps_tt computed using values of the u and d
quark parameters compatible at 2sigma with LMA+coherent from 1805.04530.
Units: [Adimensional]
"""

EPS_2 = [EPS_EE, EPS_EM, EPS_MM]
r"""float: Module-level constant

Vector of total NSI strength parameters for two-neutrino oscillations.
Units: [Adimensional]
"""

EPS_3 = [EPS_EE, EPS_EM, EPS_ET, EPS_MM, EPS_MT, EPS_TT]
r"""float: Module-level constant

Vector of total NSI strength parameters for three-neutrino oscillations.
Units: [Adimensional]
"""

# LIV parameters
# Compatible with 90% C.L. upper limits on c^(4) from 1709.03434
SXI12 = 0.0
r"""float: Module-level constant

LIV lepton mixing angle sin(xi_12).
Units: [Adimensional]
"""

SXI23 = 0.0
r"""float: Module-level constant

LIV lepton mixing angle sin(xi_23).
Units: [Adimensional]
"""

SXI13 = 0.0
r"""float: Module-level constant

LIV lepton mixing angle sin(xi_13).
Units: [Adimensional]
"""

DXICP = 0.0
r"""float: Module-level constant

LIV CP-violation phase.
Units: [radian]
"""

B1 = 1.e-9
r"""float: Module-level constant

LIV eigenvalue b_1.
Units: [eV]
"""

B2 = 1.e-9
r"""float: Module-level constant

LIV eigenvalue b_2.
Units: [eV]
"""
B3 = 2.e-9
r"""float: Module-level constant

LIV eigenvalue b_3.
Units: [eV]
"""

LAMBDA = 1.e12 # [eV]
r"""float: Module-level constant

LIV energy scale Lambda.
Units: [eV]
"""


# =============================================================================
# Historical NuFit global-fit values (v1.0 - v6.1)
# =============================================================================
#
# NUFIT_GLOBAL_FITS collects best-fit values of the standard three-flavor
# oscillation parameters from every NuFit global-fit release, from v1.0
# (2012) to v6.1 (2025), transcribed directly from the official parameter
# tables at http://www.nu-fit.org/?q=node/12 (each release's
# "vXX.tbl-parameters.pdf"). Use load_nufit_params() below to retrieve a
# specific release/ordering/category as a plain dict with the same
# {s12, s23, s13, dCP, D21, D31} keys used throughout Magnus (e.g., directly
# splattable into any osc_prob_3nu_* function).
#
# Structure: NUFIT_GLOBAL_FITS[version]['categories'][category][ordering]
# is a dict with keys 's12', 's23', 's13' (= sin(theta_ij), not sin^2),
# 'dCP' (radian), 'D21' and 'D31' (eV^2). 'ordering' is 'NO' or 'IO'. For
# IO, 'D31' is derived as Delta m^2_32 + Delta m^2_21 (Delta m^2_32 is what
# NuFit actually tabulates for IO), matching the convention already used by
# OSC_PARAMS_NU_FIT_6_0_SK_IO above.
#
# 'category' is release-specific, since NuFit has used different secondary
# splits over its history (the first-listed category is the one
# load_nufit_params() returns by default):
#   * v1.0-v1.3: 'free_fluxes_rsbl' (reactor fluxes left free, short-baseline
#     reactor data included) vs. 'huber_fluxes_no_rsbl' (fixed Huber flux
#     model, short-baseline reactor data excluded).
#   * v2.1: 'LEM' vs. 'LID', two different treatments of the reactor
#     antineutrino anomaly/flux normalization.
#   * v4.0 onward: 'with_SK' vs. 'without_SK', whether the (non-public)
#     Super-Kamiokande atmospheric chi^2 map is folded into the fit (v6.0/6.1
#     also add the latest IceCube/DeepCore chi^2 map into the 'with_SK'
#     variant; see the NuFit 6.0/6.1 papers for details).
#   * v2.0, v2.2, v3.0-v3.2: no secondary split; a single 'default' category.
#
# v1.0-v1.3 predate NuFit's separate global fits per mass ordering: only
# Delta m^2_3l was reported for each ordering hypothesis (as Delta m^2_31
# for NO and Delta m^2_32 for IO); theta12, theta23, theta13, and :math:`\delta_\text{CP}`
# were reported as a single ordering-independent fit. For these four
# releases, 'legacy_single_fit' is True, and the same theta12/theta23/
# theta13/:math:`\delta_\text{CP}` values are stored under both 'NO' and 'IO' (only 'D31'
# differs, exactly as tabulated). Where theta23 had two disconnected
# best-fit solutions ("octants"), only the global best-fit octant is
# stored (the secondary, subleading solution is not).
#
# The Bayesian-analysis variant of NuFit 2.0 (a different statistical
# methodology, not a different data category) is not included here.

NUFIT_GLOBAL_FITS = {
    'NuFIT 1.0': {
        'year': 2012,
        'legacy_single_fit': True,
        'categories': {
            'free_fluxes_rsbl': {
                'NO': {'s12': 0.5495452666, 's23': 0.6426507605, 's13': 0.1506651917, 'dCP': 5.235987756, 'D21': 7.5e-05, 'D31': 0.002473},
                'IO': {'s12': 0.5495452666, 's23': 0.6426507605, 's13': 0.1506651917, 'dCP': 5.235987756, 'D21': 7.5e-05, 'D31': -0.002352},
            },
            'huber_fluxes_no_rsbl': {
                'NO': {'s12': 0.5576737397, 's23': 0.6449806199, 's13': 0.1596871942, 'dCP': 5.2010811709, 'D21': 7.51e-05, 'D31': 0.002489},
                'IO': {'s12': 0.5576737397, 's23': 0.6449806199, 's13': 0.1596871942, 'dCP': 5.2010811709, 'D21': 7.51e-05, 'D31': -0.0023929},
            },
        },
    },
    'NuFIT 1.1': {
        'year': 2013,
        'legacy_single_fit': True,
        'categories': {
            'free_fluxes_rsbl': {
                'NO': {'s12': 0.5531726674, 's23': 0.6610597552, 's13': 0.1519868415, 'dCP': 5.9515727493, 'D21': 7.45e-05, 'D31': 0.002421},
                'IO': {'s12': 0.5531726674, 's23': 0.6610597552, 's13': 0.1519868415, 'dCP': 5.9515727493, 'D21': 7.45e-05, 'D31': -0.0023355},
            },
            'huber_fluxes_no_rsbl': {
                'NO': {'s12': 0.5594640292, 's23': 0.6603029608, 's13': 0.1587450787, 'dCP': 6.0213859194, 'D21': 7.5e-05, 'D31': 0.002429},
                'IO': {'s12': 0.5594640292, 's23': 0.6603029608, 's13': 0.1587450787, 'dCP': 6.0213859194, 'D21': 7.5e-05, 'D31': -0.002347},
            },
        },
    },
    'NuFIT 1.2': {
        'year': 2013,
        'legacy_single_fit': True,
        'categories': {
            'free_fluxes_rsbl': {
                'NO': {'s12': 0.5531726674, 's23': 0.7700649323, 's13': 0.1519868415, 'dCP': 4.6425758103, 'D21': 7.45e-05, 'D31': 0.002417},
                'IO': {'s12': 0.5531726674, 's23': 0.7700649323, 's13': 0.1519868415, 'dCP': 4.6425758103, 'D21': 7.45e-05, 'D31': -0.0023365},
            },
            'huber_fluxes_no_rsbl': {
                'NO': {'s12': 0.5594640292, 's23': 0.7694153625, 's13': 0.1562049935, 'dCP': 4.7123889804, 'D21': 7.5e-05, 'D31': 0.002429},
                'IO': {'s12': 0.5594640292, 's23': 0.7694153625, 's13': 0.1562049935, 'dCP': 4.7123889804, 'D21': 7.5e-05, 'D31': -0.002347},
            },
        },
    },
    'NuFIT 1.3': {
        'year': 2014,
        'legacy_single_fit': True,
        'categories': {
            'free_fluxes_rsbl': {
                'NO': {'s12': 0.5513619501, 's23': 0.7596051606, 's13': 0.1479864859, 'dCP': 4.3807764225, 'D21': 7.5e-05, 'D31': 0.002458},
                'IO': {'s12': 0.5513619501, 's23': 0.7596051606, 's13': 0.1479864859, 'dCP': 4.3807764225, 'D21': 7.5e-05, 'D31': -0.002373},
            },
            'huber_fluxes_no_rsbl': {
                'NO': {'s12': 0.5576737397, 's23': 0.7615773106, 's13': 0.1493318452, 'dCP': 4.5204027627, 'D21': 7.55e-05, 'D31': 0.002462},
                'IO': {'s12': 0.5576737397, 's23': 0.7615773106, 's13': 0.1493318452, 'dCP': 4.5204027627, 'D21': 7.55e-05, 'D31': -0.0023775},
            },
        },
    },
    'NuFIT 2.0': {
        'year': 2014,
        'legacy_single_fit': False,
        'categories': {
            'default': {
                'NO': {'s12': 0.5513619501, 's23': 0.6723094526, 's13': 0.1476482306, 'dCP': 5.3407075111, 'D21': 7.5e-05, 'D31': 0.002457},
                'IO': {'s12': 0.5513619501, 's23': 0.7609204952, 's13': 0.1479864859, 'dCP': 4.4331363001, 'D21': 7.5e-05, 'D31': -0.002374},
            },
        },
    },
    'NuFIT 2.1': {
        'year': 2016,
        'legacy_single_fit': False,
        'categories': {
            'LEM': {
                'NO': {'s12': 0.554977477, 's23': 0.757627877, 's13': 0.1473091986, 'dCP': 4.7472955654, 'D21': 7.49e-05, 'D31': 0.002484},
                'IO': {'s12': 0.554977477, 's23': 0.7609204952, 's13': 0.1486606875, 'dCP': 4.4680428851, 'D21': 7.49e-05, 'D31': -0.0023921},
            },
            'LID': {
                'NO': {'s12': 0.554977477, 's23': 0.6715653356, 's13': 0.1479864859, 'dCP': 5.2883476335, 'D21': 7.49e-05, 'D31': 0.002477},
                'IO': {'s12': 0.554977477, 's23': 0.7589466384, 's13': 0.1479864859, 'dCP': 4.5727626402, 'D21': 7.49e-05, 'D31': -0.0023901},
            },
        },
    },
    'NuFIT 2.2': {
        'year': 2016,
        'legacy_single_fit': False,
        'categories': {
            'default': {
                'NO': {'s12': 0.554977477, 's23': 0.6633249581, 's13': 0.1470714112, 'dCP': 5.0440015383, 'D21': 7.49e-05, 'D31': 0.002526},
                'IO': {'s12': 0.554977477, 's23': 0.764198927, 's13': 0.147478812, 'dCP': 4.6949356879, 'D21': 7.49e-05, 'D31': -0.0024431},
            },
        },
    },
    'NuFIT 3.0': {
        'year': 2016,
        'legacy_single_fit': False,
        'categories': {
            'default': {
                'NO': {'s12': 0.5531726674, 's23': 0.6640783086, 's13': 0.1471733672, 'dCP': 4.5553093477, 'D21': 7.5e-05, 'D31': 0.002524},
                'IO': {'s12': 0.5531726674, 's23': 0.7661592524, 's13': 0.1476143624, 'dCP': 4.834562028, 'D21': 7.5e-05, 'D31': -0.002439},
            },
        },
    },
    'NuFIT 3.1': {
        'year': 2017,
        'legacy_single_fit': False,
        'categories': {
            'default': {
                'NO': {'s12': 0.5540758071, 's23': 0.7516648189, 's13': 0.1481553239, 'dCP': 3.9793506945, 'D21': 7.4e-05, 'D31': 0.002515},
                'IO': {'s12': 0.5540758071, 's23': 0.756306816, 's13': 0.1487279395, 'dCP': 4.9043751981, 'D21': 7.4e-05, 'D31': -0.002409},
            },
        },
    },
    'NuFIT 3.2': {
        'year': 2018,
        'legacy_single_fit': False,
        'categories': {
            'default': {
                'NO': {'s12': 0.5540758071, 's23': 0.7334848328, 's13': 0.148526092, 'dCP': 4.0840704497, 'D21': 7.4e-05, 'D31': 0.002494},
                'IO': {'s12': 0.5540758071, 's23': 0.7443117626, 's13': 0.149231364, 'dCP': 4.8520153205, 'D21': 7.4e-05, 'D31': -0.002391},
            },
        },
    },
    'NuFIT 4.0': {
        'year': 2018,
        'legacy_single_fit': False,
        'categories': {
            'with_SK': {
                'NO': {'s12': 0.5567764363, 's23': 0.7628892449, 's13': 0.1496662955, 'dCP': 3.7873644768, 'D21': 7.39e-05, 'D31': 0.002525},
                'IO': {'s12': 0.5567764363, 's23': 0.7628892449, 's13': 0.1504327092, 'dCP': 4.8869219056, 'D21': 7.39e-05, 'D31': -0.0024381},
            },
            'without_SK': {
                'NO': {'s12': 0.5567764363, 's23': 0.7615773106, 's13': 0.1496996994, 'dCP': 3.7524578918, 'D21': 7.39e-05, 'D31': 0.002525},
                'IO': {'s12': 0.5567764363, 's23': 0.764198927, 's13': 0.150465943, 'dCP': 4.9567350757, 'D21': 7.39e-05, 'D31': -0.0024381},
            },
        },
    },
    'NuFIT 4.1': {
        'year': 2019,
        'legacy_single_fit': False,
        'categories': {
            'with_SK': {
                'NO': {'s12': 0.5567764363, 's23': 0.7503332593, 's13': 0.1495660389, 'dCP': 3.8571776469, 'D21': 7.39e-05, 'D31': 0.002528},
                'IO': {'s12': 0.5567764363, 's23': 0.7516648189, 's13': 0.1502997006, 'dCP': 4.9218284906, 'D21': 7.39e-05, 'D31': -0.0024361},
            },
            'without_SK': {
                'NO': {'s12': 0.5567764363, 's23': 0.7469939759, 's13': 0.1496996994, 'dCP': 3.8746309394, 'D21': 7.39e-05, 'D31': 0.002523},
                'IO': {'s12': 0.5567764363, 's23': 0.7503332593, 's13': 0.1503662196, 'dCP': 4.9741883682, 'D21': 7.39e-05, 'D31': -0.0024351},
            },
        },
    },
    'NuFIT 5.0': {
        'year': 2020,
        'legacy_single_fit': False,
        'categories': {
            'with_SK': {
                'NO': {'s12': 0.5513619501, 's23': 0.7569676347, 's13': 0.1489630827, 'dCP': 3.4382986264, 'D21': 7.42e-05, 'D31': 0.002517},
                'IO': {'s12': 0.5513619501, 's23': 0.7582875444, 's13': 0.1495994652, 'dCP': 4.9218284906, 'D21': 7.42e-05, 'D31': -0.0024238},
            },
            'without_SK': {
                'NO': {'s12': 0.5513619501, 's23': 0.7549834435, 's13': 0.1490301983, 'dCP': 3.4033920414, 'D21': 7.42e-05, 'D31': 0.002514},
                'IO': {'s12': 0.5513619501, 's23': 0.7582875444, 's13': 0.1496662955, 'dCP': 4.9916416607, 'D21': 7.42e-05, 'D31': -0.0024228},
            },
        },
    },
    'NuFIT 5.1': {
        'year': 2021,
        'legacy_single_fit': False,
        'categories': {
            'with_SK': {
                'NO': {'s12': 0.5513619501, 's23': 0.6708203932, 's13': 0.1498666074, 'dCP': 4.0142572796, 'D21': 7.42e-05, 'D31': 0.00251},
                'IO': {'s12': 0.5513619501, 's23': 0.7549834435, 's13': 0.1496996994, 'dCP': 4.8520153205, 'D21': 7.42e-05, 'D31': -0.0024158},
            },
            'without_SK': {
                'NO': {'s12': 0.5513619501, 's23': 0.7569676347, 's13': 0.1489966443, 'dCP': 3.3859387489, 'D21': 7.42e-05, 'D31': 0.002515},
                'IO': {'s12': 0.5513619501, 's23': 0.7602631123, 's13': 0.1495994652, 'dCP': 5.0090949532, 'D21': 7.42e-05, 'D31': -0.0024238},
            },
        },
    },
    'NuFIT 5.2': {
        'year': 2022,
        'legacy_single_fit': False,
        'categories': {
            'with_SK': {
                'NO': {'s12': 0.5504543578, 's23': 0.6715653356, 's13': 0.1491643389, 'dCP': 4.0491638646, 'D21': 7.41e-05, 'D31': 0.002507},
                'IO': {'s12': 0.5504543578, 's23': 0.7543208866, 's13': 0.1490972837, 'dCP': 4.8171087355, 'D21': 7.41e-05, 'D31': -0.0024119},
            },
            'without_SK': {
                'NO': {'s12': 0.5504543578, 's23': 0.756306816, 's13': 0.1484250653, 'dCP': 3.4382986264, 'D21': 7.41e-05, 'D31': 0.002511},
                'IO': {'s12': 0.5504543578, 's23': 0.7602631123, 's13': 0.1489630827, 'dCP': 4.9916416607, 'D21': 7.41e-05, 'D31': -0.0024239},
            },
        },
    },
    'NuFIT 5.3': {
        'year': 2024,
        'legacy_single_fit': False,
        'categories': {
            'with_SK': {
                'NO': {'s12': 0.5540758071, 's23': 0.6737952211, 's13': 0.1491308151, 'dCP': 4.0491638646, 'D21': 7.41e-05, 'D31': 0.002505},
                'IO': {'s12': 0.5540758071, 's23': 0.7536577473, 's13': 0.1490637448, 'dCP': 4.7647488579, 'D21': 7.41e-05, 'D31': -0.0024129},
            },
            'without_SK': {
                'NO': {'s12': 0.5540758071, 's23': 0.756306816, 's13': 0.1484250653, 'dCP': 3.4382986264, 'D21': 7.41e-05, 'D31': 0.002511},
                'IO': {'s12': 0.5540758071, 's23': 0.7602631123, 's13': 0.1489630827, 'dCP': 4.9916416607, 'D21': 7.41e-05, 'D31': -0.0024239},
            },
        },
    },
    'NuFIT 6.0': {
        'year': 2024,
        'legacy_single_fit': False,
        'categories': {
            'with_SK': {
                'NO': {'s12': 0.554977477, 's23': 0.68556546, 's13': 0.1488287607, 'dCP': 3.7000980142, 'D21': 7.49e-05, 'D31': 0.002513},
                'IO': {'s12': 0.554977477, 's23': 0.7416198487, 's13': 0.149365324, 'dCP': 4.7822021505, 'D21': 7.49e-05, 'D31': -0.0024091},
            },
            'without_SK': {
                'NO': {'s12': 0.5540758071, 's23': 0.7489993324, 's13': 0.1481553239, 'dCP': 3.089232776, 'D21': 7.49e-05, 'D31': 0.002534},
                'IO': {'s12': 0.5540758071, 's23': 0.7496665926, 's13': 0.1491308151, 'dCP': 4.9741883682, 'D21': 7.49e-05, 'D31': -0.0024351},
            },
        },
    },
    'NuFIT 6.1': {
        'year': 2025,
        'legacy_single_fit': False,
        'categories': {
            'with_SK': {
                'NO': {'s12': 0.5556977596, 's23': 0.68556546, 's13': 0.1499333185, 'dCP': 3.7000980142, 'D21': 7.537e-05, 'D31': 0.002511},
                'IO': {'s12': 0.5556977596, 's23': 0.7416198487, 's13': 0.1503994681, 'dCP': 4.7822021505, 'D21': 7.537e-05, 'D31': -0.00240763},
            },
            'without_SK': {
                'NO': {'s12': 0.5556977596, 's23': 0.68556546, 's13': 0.149966663, 'dCP': 3.6128315516, 'D21': 7.537e-05, 'D31': 0.002521},
                'IO': {'s12': 0.5556977596, 's23': 0.7449832213, 's13': 0.1503662196, 'dCP': 4.9392817831, 'D21': 7.537e-05, 'D31': -0.00242463},
            },
        },
    },
}


def load_nufit_params(version='NuFIT 6.1', ordering='NO', category=None):
    r"""Load standard three-flavor mixing parameters from a NuFit global fit.

    Looks up ``NUFIT_GLOBAL_FITS`` for the requested release, mass
    ordering, and (release-specific) secondary category, and returns them
    as a plain dict with the same parameter names used throughout Magnus
    (``s12``, ``s23``, ``s13``, ``dCP``, ``D21``, ``D31``), so the result
    can be passed directly as keyword arguments to any ``osc_prob_3nu_*``
    function (or to :func:`magnus.hamiltonians.hamiltonians3nu` and other
    functions that take the same standard-oscillation parameter names).

    Parameters
    ----------
    version : str, optional
        NuFit release to load, e.g. ``'NuFIT 6.1'``, ``'NuFIT 5.2'``,
        ``'NuFIT 1.0'``. See ``NUFIT_GLOBAL_FITS.keys()`` for the full list
        of available releases (v1.0 through v6.1). Default: ``'NuFIT 6.1'``
        (the latest release at the time of writing).
    ordering : str, optional
        Neutrino mass ordering: ``'NO'`` (normal) or ``'IO'`` (inverted).
        Default: ``'NO'``.
    category : str or None, optional
        Release-specific secondary category (e.g. ``'with_SK'`` /
        ``'without_SK'`` for v4.0+, ``'LEM'`` / ``'LID'`` for v2.1,
        ``'free_fluxes_rsbl'`` / ``'huber_fluxes_no_rsbl'`` for v1.0-v1.3).
        If ``None`` (default), the release's preferred/primary category is
        used (for releases with a ``with_SK``/``without_SK`` split, this is
        ``'with_SK'``). See ``NUFIT_GLOBAL_FITS[version]['categories'].keys()``
        for the categories available for a given release.

    Returns
    -------
    dict
        Dict with keys ``s12``, ``s23``, ``s13`` (:math:`\sin\theta_{ij}`,
        adimensional), ``dCP`` (radian), ``D21`` and ``D31`` (eV^2).

    Raises
    ------
    ValueError
        If ``version`` is not a known NuFit release, if ``ordering`` is
        not ``'NO'`` or ``'IO'``, or if ``category`` is not one of the
        categories available for ``version``.

    Examples
    --------
    >>> import magnus.globaldefs as gd
    >>> import magnus.oscprob as oscprob
    >>> params = gd.load_nufit_params('NuFIT 6.1', ordering='NO')
    >>> sorted(params.keys())
    ['D21', 'D31', 'dCP', 's12', 's13', 's23']
    >>> oscprob.osc_prob_3nu_vacuum(1.0 * gd.UNIT_GEV, 100.0 * gd.UNIT_KM,
    ...                             **params) # doctest: +SKIP

    Comparing normal- and inverted-ordering best fits from an older release:

    >>> no = gd.load_nufit_params('NuFIT 4.0', ordering='NO', category='with_SK')
    >>> io = gd.load_nufit_params('NuFIT 4.0', ordering='IO', category='with_SK')
    >>> no['D31'] > 0 and io['D31'] < 0
    True

    .. versionadded:: 0.11.0
    """
    if version not in NUFIT_GLOBAL_FITS:
        available = ', '.join(NUFIT_GLOBAL_FITS.keys())
        raise ValueError(
            "Error in magnus: globaldefs.load_nufit_params: unknown NuFit "
            "version '%s'. Available versions: %s." % (version, available))

    if ordering not in ('NO', 'IO'):
        raise ValueError(
            "Error in magnus: globaldefs.load_nufit_params: ordering must "
            "be 'NO' or 'IO', got '%s'." % ordering)

    categories = NUFIT_GLOBAL_FITS[version]['categories']

    if category is None:
        category = next(iter(categories))
    elif category not in categories:
        raise ValueError(
            "Error in magnus: globaldefs.load_nufit_params: category '%s' "
            "is not available for %s. Available categories: %s."
            % (category, version, ', '.join(categories.keys())))

    return dict(categories[category][ordering])


name = 'globaldefs'

__all__ = [s for s in dir() if not s.startswith('_')]
