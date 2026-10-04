# -*- coding: utf-8 -*-
r"""Notebook 28's long-range solar figure reads the probability out where a detector would.

The figure adds an :math:`L_e - L_\mu` potential, sourced by the Sun's electrons through a
light mediator, to the solar Hamiltonian.  Until this test, the sweep ended at the edge of
the solar table, 0.98 R_sun, and the averaged probability was read out in the eigenbasis of
the Hamiltonian there.  At 1/m = R_sun the new potential is still about 1e3 V_CC at that
edge, and it tilts the eigenbasis by an amount that grows with the energy: the curve turned
up above about 9 MeV and crossed the standard one.  No detector sees that.  The potential
fades over a few solar radii outside the Sun, the exit is adiabatic, and on Earth the
readout is the vacuum one.

The sweep now starts at 0.05 R_sun, where 8B neutrinos are made, ends at 20 R_sun, and
uses one coupling for both mediator ranges.  These tests lift the cell's own code out of
the generator, as `test_paper_cache_only.py` does, so they hold the notebook itself rather
than a copy of it.
"""

import pathlib
import warnings

import numpy as np
import pytest

import magnus.globaldefs as gd
import magnus.hamiltonians as hamiltonians
import magnus.matter as matter
import magnus.oscprob as oscprob
import magnus.solarmodels as solarmodels

# Reads notebooks/make_notebooks.py, which the sdist does not ship (issue #164 §1).
pytestmark = pytest.mark.checkout_only

ROOT = pathlib.Path(__file__).resolve().parent.parent
START = '# ------------------------------------------------ Figure 5b: L_e - L_mu in the Sun'
STOP = "_got = cached('solar_long_range',"


def _quiet(f, *args, **kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return f(*args, **kwargs)


@pytest.fixture(scope='module')
def cell():
    """The long-range cell of notebook 28, up to its cache call, run on the shipped model."""
    source = (ROOT/'notebooks'/'make_notebooks.py').read_text()
    assert START in source and STOP in source, (
        "make_notebooks.py no longer has notebook 28's long-range cell where this test "
        'expects it; update the anchors')
    start = source.index(START)
    code = source[start:source.index(STOP, start)]
    osc = gd.load_nufit_params('NuFIT 6.1', 'NO')
    x_solar = (solarmodels.load_solar_model('BS05-AGS-OP')['r_over_r_sun']
               * gd.SUN_RADIUS*gd.UNIT_KM)
    ne = solarmodels.electron_density_profile('BS05-AGS-OP')
    per_ne = matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
    # The notebook's names for the shipped profile, including those the cell read before
    # this fix (ne_tab, ne_sun, VCC0), so that the old cell runs and fails on its physics.
    ns = dict(np=np, gd=gd, hamiltonians=hamiltonians, oscprob=oscprob, quiet=_quiet,
              RED='red', BLUE='blue', x_solar=x_solar, R_SUN=float(x_solar[-1]),
              _ne_pkg=ne, ne_sun=ne, ne_tab=ne(x_solar), PER_NE=per_ne,
              VCC0=float(per_ne*ne(0.0)),
              HV3=np.asarray(hamiltonians.hamiltonian_3nu_vacuum_energy_independent(**osc),
                             dtype=complex))
    exec(compile(code, '<notebook 28, Figure 5b>', 'exec'), ns)
    return ns


def test_the_neutrino_is_made_where_8b_neutrinos_are(cell):
    # 8B production peaks at 0.044 R_sun in the BS05(AGS,OP) flux table; its median is 0.048.
    assert 0.04 <= cell['LR_R0']/cell['R_SUN'] <= 0.06


def test_both_ranges_share_one_coupling(cell):
    couplings = {cell['G2'][frac] for frac, _, _ in cell['LR_RANGES']}
    assert len(couplings) == 1


def test_the_readout_is_in_vacuum(cell):
    """Where the sweep ends, the new potential is negligible against the smallest vacuum
    splitting at the highest energy drawn."""
    frac = max(f for f, _, _ in cell['LR_RANGES'])
    v_end = cell['G2'][frac]*np.interp(cell['LR_END'], cell['X_LR'], cell['V_LR'][frac])
    osc = gd.load_nufit_params('NuFIT 6.1', 'NO')
    splitting = osc['D21']/(2.0*float(np.max(cell['E_LR'])))
    assert v_end < 1.0e-6*splitting


def test_the_long_range_curve_does_not_turn_up_at_high_energy(cell):
    """At 20 MeV, the longer range lowers the averaged probability below the standard one.
    Read out at the table's edge, it raised it by 0.02."""
    cell['E_LR'] = np.array([20.0])*gd.UNIT_MEV
    got = cell['_lri_sweep']()
    assert got['1'][0] < got['std'][0]
