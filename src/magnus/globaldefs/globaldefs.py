# -*- coding: utf-8 -*-
r"""Contains physical constants and unit-conversion constants.

This module contains contains values of physical constants and
unit-conversion factors used by the various modules of NuOscProbExact.
The core modules oscprob2nu.py and oscprob3nu.py do not require these
constants.

Created: 2019/04/17 17:03
Last modified: 2024/12/11 20:32
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
Units: [eV^{-1}]
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

Lepton mixing angle sin(theta_12), best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [Adimensional]
"""

S23_NO_BF_NUFIT_6_0 = np.sqrt(0.470)
r"""float: Module-level constant

Lepton mixing angle sin(theta_23), best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [Adimensional]
"""

S13_NO_BF_NUFIT_6_0 = np.sqrt(2.215e-2)
r"""float: Module-level constant

Lepton mixing angle sin(theta_13), best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [Adimensional]
"""

DCP_NO_BF_NUFIT_6_0 = 212./180.*np.pi
r"""float: Module-level constant

Lepton CP-violation phase delta_CP, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [radian]
"""

D21_NO_BF_NUFIT_6_0 = 7.49e-5
r"""float: Module-level constant

Mass-squared difference Delta m^2_21, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

D31_NO_BF_NUFIT_6_0 = 2.513e-3
r"""float: Module-level constant

Mass-squared difference Delta m^2_31, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

S12_IO_BF_NUFIT_6_0 = np.sqrt(0.308)
r"""float: Module-level constant

Lepton mixing angle sin(theta_12), best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [Adimensional]
"""

S23_IO_BF_NUFIT_6_0 = np.sqrt(0.550)
r"""float: Module-level constant

Lepton mixing angle sin(theta_23), best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [Adimensional]
"""

S13_IO_BF_NUFIT_6_0 = np.sqrt(2.231e-2)
r"""float: Module-level constant

Lepton mixing angle sin(theta_13), best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [Adimensional]
"""

DCP_IO_BF_NUFIT_6_0 = 274./180.*np.pi
r"""float: Module-level constant

Lepton CP-violation phase delta_CP, best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [radian]
"""

D21_IO_BF_NUFIT_6_0 = 7.49e-5
r"""float: Module-level constant

Mass-squared difference Delta m^2_21, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

D32_IO_BF_NUFIT_6_0 = -2.484e-3
r"""float: Module-level constant

Mass-squared difference Delta m^2_32, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [eV^2]
"""

D31_IO_BF_NUFIT_6_0 = D32_IO_BF_NUFIT_6_0+D21_IO_BF_NUFIT_6_0
r"""float: Module-level constant

Mass-squared difference Delta m^2_31, best fit from NuFit 6.0, assuming
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
Used in oscprob2nu_plot.py.
Units: [Adimensional]
"""

EPS_3 = [EPS_EE, EPS_EM, EPS_ET, EPS_MM, EPS_MT, EPS_TT]
r"""float: Module-level constant

Vector of total NSI strength parameters for three-neutrino oscillations.
Used in oscprob3nu_plot.py.
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
