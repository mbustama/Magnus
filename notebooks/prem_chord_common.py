# -*- coding: utf-8 -*-
r"""prem_chord_common.py

The PREM chord shared by the Earth analogue of Figure 11 and by its reference.

Both the reference builder and the benchmark runner have to be given *the same*
Hamiltonian, or the comparison measures bookkeeping rather than physics.  Defining the
chord once, here, is what makes that true by construction rather than by inspection.

THE CHORD.  cos(theta_z) = -0.9, so 11467.8 km of Earth, crossing sixteen PREM layer
boundaries.  Unlike the exponential profile of Figure 11, this one is *piecewise*
smooth: polynomial in radius inside each shell, discontinuous between them.  That
distinction is the whole reason this file exports the edges as well as the potential --
see ``gen_prem_reference.py`` for why a reference that ignores them is unsound.

UNITS.  ``earth.distance_traveled_inside_earth`` returns kilometres and the probability
API wants eV^-1; handing over the raw number produces a converged, unitary, wrong
answer.  The conversion happens once, here.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import numpy as np

from magnus import earth, globaldefs as gd, matter

COSTHZ = -0.9
ELECTRON_FRACTION = 0.5


def chord():
    """Returns the chord's baseline (eV^-1), its layer edges (eV^-1), and V_CC(l).

    The edges are returned in the same units as the baseline so that a caller can use
    them as slab boundaries without a second conversion -- the kind of mismatch that
    would put an edge in the wrong place and be invisible in the output.
    """
    L_km = earth.distance_traveled_inside_earth(COSTHZ)
    edges_km = np.asarray(earth.prem_layer_edges_along_chord(COSTHZ), dtype=float)
    L = L_km*gd.UNIT_KM
    edges = np.unique(np.concatenate([[0.0], edges_km*gd.UNIT_KM, [L]]))

    half = 0.5*L_km
    r_min = np.sqrt(gd.EARTH_RADIUS**2 - half**2)

    def density_along_chord(l_km):
        """PREM density in g/cm^3 at distance ``l_km`` along the chord.

        The chord's closest approach to the centre is at its midpoint, by symmetry, so
        the radial distance is r_min displaced along the chord: r^2 = r_min^2 + (l-L/2)^2.
        At l = 0 and l = L this returns the Earth's radius, which is the check below.
        """
        l = np.asarray(l_km, dtype=float)
        return earth.density_matter_func_prem(np.sqrt(r_min**2 + (l - half)**2))

    def vcc(l):
        """V_CC in eV at distance ``l`` (eV^-1) along the chord."""
        scalar = np.ndim(l) == 0
        rho = np.atleast_1d(density_along_chord(np.asarray(l, dtype=float)/gd.UNIT_KM))
        out = np.array([
            matter.VCC_func(0.0, lambda _l, _n=matter.num_density_e_func(
                0.0, lambda _x, _v=v: _v, electron_fraction=ELECTRON_FRACTION,
                density_matter_is_in_g_per_cm3=True): _n)
            for v in rho])
        return float(out[0]) if scalar else out

    return dict(baseline=L, baseline_km=L_km, edges=edges,
                vcc=vcc, density_km=density_along_chord, costhz=COSTHZ)
