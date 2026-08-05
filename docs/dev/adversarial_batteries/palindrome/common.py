"""Shared setup: a real PREM chord Hamiltonian, costhz = -0.9, 3nu, NuFIT NO best fit."""
import sys
import numpy as np

sys.path.insert(0, '/home/mbustamante/Research/magnus/src')

import magnus.magnus as mg
import magnus.earth as earth
import magnus.matter as matter
import magnus.globaldefs as gd
from magnus.hamiltonians import hamiltonians3nu as h3


COSTHZ = -0.9
ENERGY = 2.0e9          # 2 GeV in eV
DCP = gd.DCP_NO_BF_NUFIT_6_0    # 212/180*pi = 3.700...


def chord_setup(costhz=COSTHZ, energy=ENERGY, dCP=DCP):
    """Returns (L, A_func, vcc_func) for the PREM chord, in natural units."""
    L = earth.distance_traveled_inside_earth(costhz)*gd.UNIT_KM   # [eV^-1]

    h_vac_ei = h3.hamiltonian_3nu_vacuum_energy_independent(
        gd.S12_NO_BF_NUFIT_6_0, gd.S23_NO_BF_NUFIT_6_0, gd.S13_NO_BF_NUFIT_6_0,
        dCP, gd.D21_NO_BF_NUFIT_6_0, gd.D31_NO_BF_NUFIT_6_0, nubar=False)

    def rho_func(l):
        return matter.num_density_e_func(
            earth.earth_radial_distance_from_depth(costhz, l/gd.UNIT_KM),
            earth.density_matter_func_prem,
            ratio_number_neutrons_to_protons=1.0, electron_fraction=0.5,
            density_matter_is_in_g_per_cm3=True)

    vcc_func = matter.vcc_func_from_rho_func(
        rho_func, 0.0, 1.0, 0.5, False, True, True)

    proj = np.zeros((3, 3))
    proj[0][0] = 1.0

    def H(l):
        vcc = np.asarray(vcc_func(l))
        return (1.0/energy)*h_vac_ei + vcc[..., None, None]*proj

    def A(l):
        return -1j*H(l)

    return L, A, H, vcc_func


def uniform_edges(L, n_slabs):
    e = np.linspace(0.0, L, n_slabs + 1)
    return np.stack([e[:-1], e[1:]], axis=1)


def symmetrised_edges(L, n_slabs):
    """Edges whose widths are exactly palindromic (w = (w + w[::-1])/2)."""
    e = np.linspace(0.0, L, n_slabs + 1)
    w = np.diff(e)
    w = 0.5*(w + w[::-1])
    a = np.concatenate([[0.0], np.cumsum(w)])
    return np.stack([a[:-1], a[1:]], axis=1), w


def sample_A(A, edges, n_tpts, method, order):
    """Samples of A on each slab, shape (n_slabs, m, d, d), plus widths."""
    edges = np.asarray(edges, dtype=float)
    widths = edges[:, 1] - edges[:, 0]
    if method == 'gl':
        s = mg.gl_nodes(order)
    else:
        s = np.linspace(0.0, 1.0, n_tpts)
    tgrid = edges[:, :1] + widths[:, None]*s
    At = A(tgrid.ravel()).reshape(edges.shape[0], len(s), 3, 3)
    return At, widths
