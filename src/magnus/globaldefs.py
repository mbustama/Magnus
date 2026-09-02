# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""globaldefs.py

Contains physical constants and unit-conversion constants.

This module contains values of physical constants and unit-conversion
factors used by the various modules of Magnus: unit conversions (km,
cm, GeV, etc., to natural units of eV), fundamental constants (G_F,
particle masses, Avogadro's number), Earth/Sun radii and reference
densities, flavor index constants (NUE, NUMU, NUTAU, NUS), predefined
oscillation/NSI/LIV parameter sets (e.g., NuFit 6.1, the default), and ANSI terminal
color codes (class ``cstyle``) used to format warning/error messages.

Routine listings
----------------

    * cstyle - ANSI terminal color-code constants
    * set_color_output - Enables or disables ANSI color in the warning
           and error message prefixes
    * load_nufit_params - Loads one NuFit release/ordering/category as a
           dict of standard oscillation parameters, in whichever ``angles``
           convention is asked for

The remaining module-level names are physical constants, unit-conversion
factors, the ANGLE_CONVENTIONS tuple and the MixingAngleConventionWarning
class, not routines; see the module source for the full list.
"""


__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


# from numpy import *
import numpy as np

import os
import platform

# If on Windows, need to call os.system() to print in color in stdout
if platform.system() == 'Windows':  # pragma: no cover - Windows only
    os.system("")

# Class of different styles
class cstyle():
    r"""ANSI escape-code constants for colored/styled terminal output.

    Used to format the warning/error/tolerance messages printed by
    ``oscprob.py`` (e.g., ``gd.WARNING_MSG_IN_COLOR``,
    ``gd.ERROR_MSG_IN_COLOR``). Has no effect on Windows terminals unless
    ``os.system("")`` has been called first, which this module does at
    import time.

    .. versionadded:: 1.0.0
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

ANGLE_CONVENTIONS = ('sin', 'sin2', 'rad', 'deg')
r"""tuple: The values the ``angles`` keyword accepts, in the order they are documented.

``'sin'`` (the default everywhere) is the sine of the mixing angle, ``'sin2'`` its
square -- which is what global fits report -- ``'rad'`` the angle itself in radians, and
``'deg'`` in degrees.  Under ``'deg'`` the CP phases are read as degrees too; under the
other three they stay in radians, a sine being no way to state a phase.

Defined here rather than in :mod:`magnus.hamiltonians` because the conversion itself lives
in a private module, and a name users are told to filter or compare against has to be
documented somewhere public.

.. versionadded:: 1.0.0
"""


class MixingAngleConventionWarning(UserWarning):
    r"""A parameter set is very probably not in the ``angles`` convention it declared.

    Raised only where the mistake is diagnosable from the values themselves: sines handed
    to ``angles='deg'`` are about fifty times too small to be angles, and the call would
    otherwise return a converged, unitary, entirely wrong probability rather than an
    error.  A warning rather than an exception, for the same reason as
    :class:`magnus.matter.DensityUnitWarning`: the threshold reflects the mixing people
    currently study, not a law.

    Its own class so it can be silenced or promoted on its own::

        import warnings
        import magnus.globaldefs as gd

        warnings.filterwarnings('error', category=gd.MixingAngleConventionWarning)

    .. versionadded:: 1.0.0
    """


class BaselineUnitWarning(UserWarning):
    r"""A baseline was passed that looks like kilometers rather than eV\ :sup:`-1`.

    Every length crossing this API is in natural units, so a baseline is
    :math:`L_{\rm km} \times` :data:`CONV_KM_TO_INV_EV`, some 5.07e9 per kilometer.
    Passing the raw kilometer value does not fail: the call returns a converged, exactly
    unitary probability for a baseline a few meters long, which looks like an ordinary
    answer rather than a wrong one.  Measured on the Sun, 694700 passed raw returns 0.910
    at 20 MeV where the correct value is 0.290, and the survival probability comes out
    *rising* with energy, which is backwards for an MSW resonance.

    The threshold is :data:`IMPLAUSIBLE_BASELINE_NATURAL_UNITS`, about two meters in
    natural units, so a genuinely short baseline is still reachable without tripping it.

    Its own class so it can be silenced deliberately::

        import warnings
        import magnus.globaldefs as gd

        warnings.filterwarnings('ignore', category=gd.BaselineUnitWarning)

    .. versionadded:: 1.0.5
    """


#: Below this, a baseline in natural units is almost certainly kilometers left
#: unconverted: 1e7 eV^-1 is about two meters, and the shortest oscillation baselines
#: anyone runs are meters, which land above it once converted.
IMPLAUSIBLE_BASELINE_NATURAL_UNITS = 1.0e7


class SterileMatterCompositionWarning(UserWarning):
    r"""The sterile matter entry is built from a different medium than the density.

    An Earth chord takes its neutron-to-proton ratio from :math:`Y_e` layer by layer for the
    density, but the sterile states' entry in the matter projector is one matrix for the
    whole chord and takes the caller's scalar instead.  They disagree by construction unless
    the caller matches them, and the disagreement is worth about 2e-02 in probability at 3+1
    on a core-crossing chord -- silently, since nothing else about the call looks wrong.

    Its own class so it can be silenced once the choice has been made deliberately::

        import warnings
        import magnus.globaldefs as gd

        warnings.filterwarnings('ignore', category=gd.SterileMatterCompositionWarning)

    Three flavors never raise it: the projector's sterile block is empty.

    .. versionadded:: 1.0.0
    """


TOL_MSG_NO_COLOR = "Requested tolerance achieved"

TOL_MSG_IN_COLOR = cstyle.CGREENBG + "Requested tolerance achieved" + cstyle.CEND


def set_color_output(enabled: bool) -> None:
    r"""Enables or disables ANSI color in the warning/error/tolerance message prefixes.

    The ``*_IN_COLOR`` constants above wrap their text in ANSI escape codes, which read
    correctly in a terminal but appear as literal escape-code noise anywhere that does not
    interpret them -- a captured log file, a Jupyter notebook rendered to HTML, or the
    ``jupyter-execute`` cells in this package's own documentation.  Calling this function with
    ``False`` rebinds every ``*_IN_COLOR`` constant to its plain-text counterpart, so the
    existing call sites (which all reference the ``*_IN_COLOR`` names) print unadorned text
    with no further change.  Calling it with ``True`` restores the colored versions.

    .. versionadded:: 1.0.0

    Parameters
    ----------
    enabled : bool
        True to emit ANSI-colored message prefixes (the default at import time), False to emit
        plain text.

    Returns
    -------
    None

    Examples
    --------
    .. jupyter-execute::

        import magnus.globaldefs as gd

        gd.set_color_output(False)
        gd.WARNING_MSG_IN_COLOR
    """
    global WARNING_MSG_IN_COLOR, ERROR_MSG_IN_COLOR, TOL_MSG_IN_COLOR
    if enabled:
        WARNING_MSG_IN_COLOR = cstyle.CVIOLETBG + "Warning:" + cstyle.CEND
        ERROR_MSG_IN_COLOR = cstyle.CREDBG + "Error in magnus:" + cstyle.CEND
        TOL_MSG_IN_COLOR = cstyle.CGREENBG + "Requested tolerance achieved" + cstyle.CEND
    else:
        WARNING_MSG_IN_COLOR = WARNING_MSG_NO_COLOR
        ERROR_MSG_IN_COLOR = ERROR_MSG_NO_COLOR
        TOL_MSG_IN_COLOR = TOL_MSG_NO_COLOR


MAGNUS_MAX_PREDEFINED_NUM_FLAVORS = 5
r"""float: Module-level constant

Maximum number of flavors for which we have hard-coded routines in the oscprob module.
Units: [Adimensional]
"""


from magnus.magnus import MAGNUS_EXP_ORDER_MAX  # noqa: E402,F401
r"""int: Module-level constant

Maximum order of the Magnus expansion currently supported.  Re-exported from
:mod:`magnus.magnus`, which is where the expansion is implemented and therefore the only
place the ceiling is defined; it used to be written out here as well, and the two copies
had to be kept in step by hand.
Units: [Adimensional]
"""


CONV_KM_TO_INV_EV = 5.06773e9
r"""float: Module-level constant

Multiplicative conversion factor from km to :math:`\text{eV}^{-1}`.
Units: [:math:`\text{km}^{-1}~\text{eV}^{-1}`].
"""

UNIT_KM = CONV_KM_TO_INV_EV
r"""float: Module-level constant

Multiplicative conversion factor from km to :math:`\text{eV}^{-1}`.  Alias for CONV_KM_TO_INV_EV.
Units: [:math:`\text{km}^{-1}~\text{eV}^{-1}`].
"""

CONV_CM_TO_INV_EV = CONV_KM_TO_INV_EV*1.e-5
r"""float: Module-level constant

Multiplicative conversion factor from cm to :math:`\text{eV}^{-1}`.
Units: [:math:`\text{cm}^{-1}~\text{eV}^{-1}`]
"""

UNIT_CM = CONV_CM_TO_INV_EV
r"""float: Module-level constant

Multiplicative conversion factor from cm to :math:`\text{eV}^{-1}`.  Alias for CONV_CM_TO_INV_EV.
Units: [:math:`\text{cm}^{-1}~\text{eV}^{-1}`]
"""

CONV_CM3_TO_INV_EV3 = np.power(CONV_CM_TO_INV_EV, 3.0)
r"""float: Module-level constant

Multiplicative conversion factor from :math:`\text{cm}^{3}` to :math:`\text{eV}^{-3}`.
Units: [:math:`\text{cm}^{-3}~\text{eV}^{-3}`]
"""

UNIT_CM3 = CONV_CM3_TO_INV_EV3
r"""float: Module-level constant

Multiplicative conversion factor from :math:`\text{cm}^{3}` to :math:`\text{eV}^{-3}`.  Alias for CONV_CM3_TO_INV_EV3.
Units: [:math:`\text{cm}^{-3}~\text{eV}^{-3}`]
"""

CONV_INV_EV_TO_CM = 1./CONV_CM_TO_INV_EV
r"""float: Module-level constant

Multiplicative conversion factor from :math:`\text{eV}^{-1}` to cm.
Units: [eV cm]
"""

UNIT_PER_CM3 = 1.0/CONV_CM3_TO_INV_EV3
r"""float: Module-level constant

Multiplicative conversion factor from :math:`\text{cm}^{-3}` to :math:`\text{eV}^{3}`. 
Units: [:math:`\text{cm}^{3}~\text{eV}^{3}`]
"""

CONV_EV_TO_G = 1.783e-33
r"""float: Module-level constant

Multiplicative conversion factor from :math:`\text{eV}^{-1}` to grams.
Units: [:math:`\text{g eV}^{-1}`]
"""

CONV_G_TO_EV = 1./CONV_EV_TO_G
r"""float: Module-level constant

Multiplicative conversion factor from grams to :math:`\text{eV}^{-1}`.
Units: [:math:`\text{eV g}^{-1}`]
"""

UNIT_G_PER_CM3 = CONV_G_TO_EV/CONV_CM3_TO_INV_EV3
r"""float: Module-level constant

Multiplicative conversion factor from :math:`\text{g cm}^{-3}` to :math:`\text{eV}^{4}`.
Units: [:math:`\text{g}^{-1}~\text{cm}^{3}~\text{eV}^{4}`]
"""


SQRT_OF_2 = np.sqrt(2.0)
r"""float: Module-level constant

Square root of 2..
Units: [Adimensional]
"""

GF = 1.1663787e-23
r"""float: Module-level constant

Fermi constant.
Units: [:math:`\text{eV}^{-2}`]
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
Units: [:math:`\text{g cm}^{-3}`]
"""

N_AV = 6.02214076e23
r"""float: Module-level constant

Avogadro constant
Units: [:math:`\text{mol}^{-1}`]
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
Units: [:math:`\text{eV}^{3}`]
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

Normalization of the standard exponential fit to the solar electron number density,
:math:`n_e(r) = 245\,N_A\,\exp(-10.54\,r/R_\odot)\ \text{cm}^{-3}`.

This is the :math:`r \to 0` **intercept of that fit**, not the central density of a solar model.
The two differ by more than a factor of two: the BS2005-AGS,OP table gives
:math:`n_e = 102.7\,N_A\ \text{cm}^{-3}` at its innermost point :math:`r = 0.0016\,R_\odot`,
against the fit's 245, because the real profile flattens towards the center while an exponential
does not.  Measured against that table across the whole star:

===========================  =========================
:math:`r/R_\odot`            largest departure from it
===========================  =========================
0.00 -- 0.05                 fit high by 2.4x
0.05 -- 0.10                 fit high by 1.6x
0.10 -- 0.20                 21 %
**0.20 -- 0.30**             **2.5 %**
0.30 -- 0.70                 11 %
0.70 -- 1.00                 up to 89 %
===========================  =========================

So the fit is a few-percent description only in a band around :math:`0.2\,R_\odot`, and the
scale height :data:`L_SCALE_SUN` -- the span the package's own diagnostics use as a trajectory --
lies at :math:`0.095\,R_\odot`, where the fit is high by about 30 %.  It is still the right
constant for the exponential profile it normalizes; what would be wrong is reading it as a
measurement of the Sun's central density, or the exponential as a stand-in for a tabulated model
in the core.

Units: [:math:`\text{eV}^{3}`]
"""

L_SCALE_SUN = SUN_RADIUS/10.54*UNIT_KM
r"""float: Module-level constant

Scale height of the solar electron number density, :math:`R_\odot/10.54`, the decay length of
the exponential fit normalized by :data:`NUM_DENSITY_E_SUN_CENTRAL`.

Units: [:math:`\text{eV}^{-1}`]
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
Units: [:math:`\text{eV}^{2}`]
"""

D31_NO_BF_NUFIT_6_0 = 2.513e-3
r"""float: Module-level constant

Mass-squared difference :math:`\Delta m_{31}^2`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [:math:`\text{eV}^{2}`]
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
Units: [:math:`\text{eV}^{2}`]
"""

D32_IO_BF_NUFIT_6_0 = -2.484e-3
r"""float: Module-level constant

Mass-squared difference :math:`\Delta m_{32}^2`, best fit from NuFit 6.0, assuming
normal ordering with SK atmospheric data.
Units: [:math:`\text{eV}^{2}`]
"""

D31_IO_BF_NUFIT_6_0 = D32_IO_BF_NUFIT_6_0+D21_IO_BF_NUFIT_6_0
r"""float: Module-level constant

Mass-squared difference :math:`\Delta m_{31}^2`, best fit from NuFit 6.0, assuming
inverted ordering with SK atmospheric data.
Units: [:math:`\text{eV}^{2}`]
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


def load_nufit_params(version='NuFIT 6.1', ordering='NO', category=None, angles='sin'):
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
    angles : str, optional
        Convention the three mixing angles are returned in: ``'sin'`` (default) their sines,
        ``'sin2'`` their sines *squared* -- which is the form NuFit itself reports --
        ``'rad'`` the angles in radians, or ``'deg'`` in degrees.  Under ``'deg'`` ``dCP``
        is converted too.

        **Pass the same value here that you pass to the probability function.**  The two are
        one setting in two places: ``osc_prob_3nu_earth(E, **load_nufit_params(), angles='deg')``
        reads perfectly and is silently wrong, because the loader's sines (0.15 to 0.85) are
        then interpreted as degrees, about fifty times too small.  The result is a converged,
        unitary, entirely wrong probability.  The guard in
        :func:`magnus.hamiltonians.hamiltonians3nu.hamiltonian_3nu_vacuum_energy_independent`
        catches that particular pairing, but the reliable fix is to state the convention once
        and use it on both calls.
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
        Dict with keys ``s12``, ``s23``, ``s13``, ``dCP``, ``D21`` and ``D31``.  The three
        angles are in whichever convention ``angles`` names -- by default their sines,
        adimensional -- and ``dCP`` is in radians unless ``angles='deg'``, which puts it in
        degrees.  ``D21`` and ``D31`` are always :math:`\text{eV}^{2}`.

    Raises
    ------
    ValueError
        If ``version`` is not a known NuFit release, if ``ordering`` is
        not ``'NO'`` or ``'IO'``, or if ``category`` is not one of the
        categories available for ``version``.

    Examples
    --------
    Load a release and feed it straight into a probability function (the
    code below runs when these docs are built, so the output shown is
    always current):

    .. jupyter-execute::

        import magnus.globaldefs as gd
        import magnus.oscprob as oscprob

        params = gd.load_nufit_params('NuFIT 6.1', ordering='NO')

        sorted(params.keys())

    .. jupyter-execute::

        oscprob.osc_prob_3nu_vacuum(1.0 * gd.UNIT_GEV, 100.0 * gd.UNIT_KM,
                                    **params)

    Comparing normal- and inverted-ordering best fits from an older release:

    .. jupyter-execute::

        no = gd.load_nufit_params('NuFIT 4.0', ordering='NO', category='with_SK')
        io = gd.load_nufit_params('NuFIT 4.0', ordering='IO', category='with_SK')

        no['D31'] > 0 and io['D31'] < 0

    .. versionadded:: 1.0.0
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

    out = dict(categories[category][ordering])

    # The stored tables are sines, so 'sin' is a pass-through and costs nothing.  The other
    # three are computed from them here rather than by the caller, because a caller who
    # converts by hand and then names a convention has two places to get it wrong.
    if angles != 'sin':
        from magnus.hamiltonians import _angles as _a
        _a.validate_convention('globaldefs.load_nufit_params', angles)
        ang, ph = _a.from_sines(angles,
                                {k: out[k] for k in ('s12', 's23', 's13')},
                                {'dCP': out['dCP']})
        out.update(ang)
        out.update(ph)

    return out


# ---------------------------------------------------------------------------------------
# The fallback set, defined here rather than beside OSC_PARAMS_PREDEFINED above because it
# is *derived* from `load_nufit_params`, which is defined immediately above this.
#
# It used to be a second copy of the numbers, and the two copies disagreed: omitting
# oscillation parameters fell back to NuFIT 6.0, while `load_nufit_params()` with no
# arguments returned 6.1.  Both were documented, so neither read as a mistake -- but the
# same script got different answers depending on which door it came through, by 4.0e-03 in
# probability at 1 GeV over 1300 km.  Deriving the fallback from the loader means there is
# one set of numbers, and a future release is a one-line change here rather than a second
# table to keep in step.
OSC_PARAMS_NU_FIT_6_1_SK_NO = {
    'name': 'OSC_PARAMS_NU_FIT_6_1_NO',
    'description': 'NuFit 6.1, NO, with SK atmospheric data',
    **load_nufit_params('NuFIT 6.1', 'NO', category='with_SK'),
}

OSC_PARAMS_NU_FIT_6_1_SK_IO = {
    'name': 'OSC_PARAMS_NU_FIT_6_1_IO',
    'description': 'NuFit 6.1, IO, with SK atmospheric data',
    **load_nufit_params('NuFIT 6.1', 'IO', category='with_SK'),
}

OSC_PARAMS_PREDEFINED['OSC_PARAMS_NU_FIT_6_1_SK_NO'] = OSC_PARAMS_NU_FIT_6_1_SK_NO
OSC_PARAMS_PREDEFINED['OSC_PARAMS_NU_FIT_6_1_SK_IO'] = OSC_PARAMS_NU_FIT_6_1_SK_IO

# The 6.0 sets stay reachable by name: this changes which release is the *default*, not
# which releases exist, so a caller pinned to 6.0 asks for it explicitly and keeps working.
OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT'] = OSC_PARAMS_NU_FIT_6_1_SK_NO


__all__ = [
    'cstyle',
    'set_color_output',
    'WARNING_MSG_NO_COLOR',
    'WARNING_MSG_IN_COLOR',
    'ERROR_MSG_NO_COLOR',
    'ERROR_MSG_IN_COLOR',
    'ANGLE_CONVENTIONS',
    'MixingAngleConventionWarning',
    'SterileMatterCompositionWarning',
    'TOL_MSG_NO_COLOR',
    'TOL_MSG_IN_COLOR',
    'MAGNUS_MAX_PREDEFINED_NUM_FLAVORS',
    'MAGNUS_EXP_ORDER_MAX',
    'CONV_KM_TO_INV_EV',
    'UNIT_KM',
    'CONV_CM_TO_INV_EV',
    'UNIT_CM',
    'CONV_CM3_TO_INV_EV3',
    'UNIT_CM3',
    'CONV_INV_EV_TO_CM',
    'UNIT_PER_CM3',
    'CONV_EV_TO_G',
    'CONV_G_TO_EV',
    'UNIT_G_PER_CM3',
    'SQRT_OF_2',
    'GF',
    'MASS_ELECTRON',
    'MASS_PROTON',
    'MASS_NEUTRON',
    'ELECTRON_FRACTION_EARTH_CRUST',
    'DENSITY_MATTER_CRUST_G_PER_CM3',
    'N_AV',
    'NUM_DENSITY_E_EARTH_CRUST',
    'VCC_EARTH_CRUST',
    'EARTH_RADIUS',
    'SUN_RADIUS',
    'NUM_DENSITY_E_SUN_CENTRAL',
    'L_SCALE_SUN',
    'NUE',
    'NUMU',
    'NUTAU',
    'NUS',
    'NUS1',
    'NUS2',
    'UNIT_KEV',
    'UNIT_MEV',
    'UNIT_GEV',
    'UNIT_TEV',
    'UNIT_PEV',
    'UNIT_EEV',
    'S12_NO_BF_NUFIT_6_0',
    'S23_NO_BF_NUFIT_6_0',
    'S13_NO_BF_NUFIT_6_0',
    'DCP_NO_BF_NUFIT_6_0',
    'D21_NO_BF_NUFIT_6_0',
    'D31_NO_BF_NUFIT_6_0',
    'S12_IO_BF_NUFIT_6_0',
    'S23_IO_BF_NUFIT_6_0',
    'S13_IO_BF_NUFIT_6_0',
    'DCP_IO_BF_NUFIT_6_0',
    'D21_IO_BF_NUFIT_6_0',
    'D32_IO_BF_NUFIT_6_0',
    'D31_IO_BF_NUFIT_6_0',
    'OSC_PARAMS_NU_FIT_6_0_SK_NO',
    'OSC_PARAMS_NU_FIT_6_0_SK_IO',
    'OSC_PARAMS_NU_FIT_6_1_SK_NO',
    'OSC_PARAMS_NU_FIT_6_1_SK_IO',
    'OSC_PARAMS_PREDEFINED',
    'EPS_EE',
    'EPS_EM',
    'EPS_ET',
    'EPS_MM',
    'EPS_MT',
    'EPS_TT',
    'EPS_2',
    'EPS_3',
    'SXI12',
    'SXI23',
    'SXI13',
    'DXICP',
    'B1',
    'B2',
    'B3',
    'LAMBDA',
    'NUFIT_GLOBAL_FITS',
    'load_nufit_params',
]
