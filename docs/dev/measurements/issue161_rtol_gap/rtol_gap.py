"""True error against the requested tolerance on smooth PREM chords (issue #161 comment).

The oracle is the one diagnostics.rst prescribes: scipy's DOP853 at rtol=1e-12, run layer by
layer between the PREM crossings, with its own movement to rtol=1e-13 printed.  Each Magnus
call is `osc_prob_3nu_earth` at 3 GeV, electron_fraction=0.5, strategy='magnus', over the
tester's grid: three quadrature methods x magnus_exp_order 2, 4, 6, 8 x rtol = atol = 1e-4
and 1e-8, on the chords at costhz = -0.8 and -1.

Run from the repository root:  python docs/dev/measurements/issue161_rtol_gap/rtol_gap.py
"""
import warnings

import numpy as np
from scipy.integrate import solve_ivp

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.hamiltonians as hh
import magnus.matter as matter
import magnus.oscprob as o

warnings.simplefilter('ignore')
p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
osc = {k: p[k] for k in ('s12', 's23', 's13', 'dCP', 'D21', 'D31')}
E = 3*gd.UNIT_GEV
hv = np.asarray(hh.hamiltonian_3nu_vacuum_energy_independent(**osc))/E
e00 = np.diag([1.0, 0.0, 0.0])


def reference(cz, tol):
    """|S|^2 from DOP853, integrated between the PREM crossings of the chord."""
    L = earth.distance_traveled_inside_earth(cz)*gd.UNIT_KM
    rho, _ = o._earth_composition(cz, 0.5, None, None, None, None, None, 'ref', num_flavors=3)
    vcc = matter.vcc_func_from_rho_func(rho, density_is_of_number_of_electrons=True)
    edges = np.concatenate([[0.0], earth.prem_layer_edges_along_chord(cz)*gd.UNIT_KM, [L]])
    edges = np.unique(edges[(edges >= 0.0) & (edges <= L)])
    psi = np.eye(3, dtype=complex)
    for a, b in zip(edges[:-1], edges[1:]):
        def rhs(l, y):
            return (-1j*(hv + float(vcc(l))*e00) @ y.reshape(3, 3)).ravel()
        sol = solve_ivp(rhs, (a, b), psi.ravel(), method='DOP853', rtol=tol, atol=1e-2*tol)
        psi = sol.y[:, -1].reshape(3, 3)
    return np.abs(psi)**2, L


rows = []
for cz in (-0.8, -1.0):
    R12, L = reference(cz, 1e-12)
    R, _ = reference(cz, 1e-13)
    print('costhz %.1f: oracle movement 1e-12 -> 1e-13 = %.1e' % (cz, np.abs(R12 - R).max()))
    for meth in ('gl', 'trapezoid', 'simpson'):
        for order in (2, 4, 6, 8):
            for tol in (1e-4, 1e-8):
                P = np.asarray(o.osc_prob_3nu_earth(
                    E, costhz=cz, L=L, electron_fraction=0.5, integration_method=meth,
                    magnus_exp_order=order, rtol=tol, atol=tol, strategy='magnus'))
                # The oracle's index order is (final, initial); compare in either layout.
                err = min(np.abs(P - R.T).max(), np.abs(P - R).max())
                rows.append((err/tol, cz, meth, order, tol, err))
rows.sort(reverse=True)
print('%d calls; above the tolerance: %d' % (len(rows), sum(r[0] > 1 for r in rows)))
print('true error / tolerance: max %.3g, median %.3g, min %.3g'
      % (rows[0][0], rows[len(rows)//2][0], rows[-1][0]))
for r in rows[:5]:
    print('  %.2f  costhz %.1f  %-9s order %d  tol %.0e  error %.3g' % r)
