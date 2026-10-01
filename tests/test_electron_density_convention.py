# -*- coding: utf-8 -*-
"""A mass density becomes an electron density as rho N_A Y_e (issue #168).

Y_e is the number of electrons per atomic mass unit of the material, sum w_i Z_i/A_i over its
elements with their atomic weights: iron's 26/55.845 and seawater's 10/18.015 are the Earth's
own core and ocean values.  With Y_e so defined, n_e = rho Y_e / m_u = rho N_A Y_e exactly.
Until 1.2.0 the package divided by the mean free-nucleon mass, (m_p + r m_n)/(1 + r), which
ignores nuclear binding and put every electron density 0.8% low, and it converted grams with a
rounded e/c^2, another 1.9e-4.  Every other code the paper benchmarks against (nuSQuIDS,
GLoBES, Prob3++, NuFast, OscProb, nuCraft) uses rho N_A Y_e.
"""

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.matter as matter
import magnus.solarmodels as solarmodels

PER_CM3 = gd.UNIT_PER_CM3


def test_the_gram_is_converted_with_the_exact_electronvolt():
    assert gd.CONV_EV_TO_G == pytest.approx(1.602176634e-19/299792458.0**2*1.0e3, rel=1e-15)


def test_one_gram_holds_avogadro_atomic_mass_units():
    # 1 g / m_u = N_A up to the 2019 SI's 3.5e-10 difference between the two.
    assert gd.CONV_G_TO_EV/gd.ATOMIC_MASS_UNIT/gd.N_AV == pytest.approx(1.0, abs=1e-9)


@pytest.mark.parametrize('ratio', [0.0, 1.0, 1.15, 2.0])
@pytest.mark.parametrize('ye', [0.4656, 0.4957, 0.5, 0.5551])
def test_the_electron_density_is_rho_avogadro_ye(ye, ratio):
    """Whatever the neutron-to-proton ratio: it no longer enters the conversion."""
    rho = 5.0
    ne = matter.num_density_e_func(0.0, lambda l: rho, ratio_number_neutrons_to_protons=ratio,
                                   electron_fraction=ye, density_matter_is_in_g_per_cm3=True)
    assert ne/(rho*gd.N_AV*ye*PER_CM3) == pytest.approx(1.0, abs=1e-9)


def test_the_potential_is_the_textbook_one():
    """V_CC = sqrt(2) G_F N_A rho Y_e = 7.6325e-14 eV per g/cm^3 of rho Y_e."""
    vcc = matter.vcc_func_from_rho_func(1.0, 0.0, 1.0, 1.0, False, True)
    assert vcc == pytest.approx(7.6325e-14, rel=1e-4)
    assert gd.VCC_EARTH_CRUST == pytest.approx(
        np.sqrt(2.0)*gd.GF*gd.DENSITY_MATTER_CRUST_G_PER_CM3*gd.ELECTRON_FRACTION_EARTH_CRUST
        *gd.N_AV*PER_CM3, rel=1e-9)


@pytest.mark.parametrize('name', ['B16-GS98', 'BS05-OP', 'BP04', 'B23-GS98'])
def test_the_solar_electron_density_counts_hydrogen_and_helium_by_their_atomic_masses(name):
    model = solarmodels.load_solar_model(name)
    rho, X = model['rho_g_per_cm3'], model['x_hydrogen']
    expected = rho*gd.N_AV*PER_CM3*(X/1.00782503 + 2.0*(1.0 - X)/4.00260325)
    r = model['r_over_r_sun']*gd.SUN_RADIUS*gd.UNIT_KM
    ne = solarmodels.electron_density_profile(name)(r)
    assert np.max(np.abs(ne/expected - 1.0)) < 1e-9
    # And the sterile states' ratio from the same composition: helium's neutrons over all protons.
    helium = 2.0*(1.0 - X)/4.00260325
    ratio = solarmodels.neutron_to_proton_ratio_profile(name)(r)
    assert np.max(np.abs(ratio - helium/(X/1.00782503 + helium))) < 1e-12


def test_the_solar_centre_moves_by_what_the_atomic_masses_say():
    """The centre of B16-GS98 (X = 0.3466): 0.37% above the mean free-nucleon count, and 0.43%
    below the textbook rho N_A (1 + X)/2 that Bahcall's electron-density tables follow."""
    model = solarmodels.load_solar_model('B16-GS98')
    X, rho = model['x_hydrogen'][0], model['rho_g_per_cm3'][0]
    ne = float(solarmodels.electron_density_profile('B16-GS98')(0.0))/(gd.N_AV*PER_CM3)
    free_nucleon = rho*(1.0 + X)/2.0/(0.5*(gd.MASS_PROTON + gd.MASS_NEUTRON)/gd.ATOMIC_MASS_UNIT)
    textbook = rho*(1.0 + X)/2.0
    assert 3e-3 < ne/free_nucleon - 1.0 < 5e-3
    assert -6e-3 < ne/textbook - 1.0 < -3e-3
