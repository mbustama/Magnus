# -*- coding: utf-8 -*-
"""Shared harness for the adversarial-validation batteries.

Every trap listed in HANDOVER_ADVERSARIAL_VALIDATION.md is encoded here once, so the
batteries cannot re-hit them:

  * vcc_func_from_rho_func's 7th positional argument is `is_number_density`, not `nubar`.
  * H closures take exactly ONE parameter (factory-built, never `lambda l, E=E: ...`).
  * solve_ivp always gets t_eval=None but a bounded max_step is NOT set; instead we solve
    for U only at the endpoint(s) we need and never store every accepted step.
  * osc_params dicts are filtered to the mixing parameters only (no name/description).
"""

import warnings

import numpy as np
from scipy.integrate import solve_ivp

import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.matter as matter

# ----------------------------------------------------------------------
# Oscillation parameters, filtered (trap: name/description are rejected downstream)
# ----------------------------------------------------------------------

_MIX_KEYS = ('s12', 's23', 's13', 'dCP', 'D21', 'D31')


def std_params():
    p = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
    return {k: p[k] for k in _MIX_KEYS}


def params_for(d, sterile_scale=1.0, seed=None):
    """Oscillation parameters for d = 2..5.  Sterile mixings are deliberately large enough
    to put extra resonances in the solar density range."""
    p = std_params()
    if d == 2:
        return {'sth': p['s12'], 'Dm2': p['D21']}
    if d == 3:
        return dict(p)
    if d == 4:
        q = dict(p)
        q.update(s14=0.30*sterile_scale, d14=0.7, s24=0.20*sterile_scale, d24=1.9,
                 s34=0.15*sterile_scale, D41=1.0e-3*sterile_scale)
        return q
    if d == 5:
        q = dict(p)
        q.update(s14=0.30*sterile_scale, d14=0.7, s15=0.25*sterile_scale, d15=2.2,
                 s24=0.20*sterile_scale, d24=1.9, s25=0.18*sterile_scale,
                 s34=0.15*sterile_scale, s35=0.12*sterile_scale, d35=0.4,
                 D41=1.0e-3*sterile_scale, D51=3.0e-3*sterile_scale)
        return q
    raise ValueError(d)


def h_vac_for(d, osc_params, nubar=False):
    p = osc_params
    if d == 2:
        return hams.hamiltonian_2nu_vacuum_energy_independent(p['sth'], p['Dm2'])
    if d == 3:
        return hams.hamiltonian_3nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['D21'], p['D31'], nubar=nubar)
    if d == 4:
        return hams.hamiltonian_4nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['s14'], p['d14'], p['s24'], p['d24'],
            p['s34'], p['D21'], p['D31'], p['D41'], nubar=nubar)
    if d == 5:
        return hams.hamiltonian_5nu_vacuum_energy_independent(
            p['s12'], p['s23'], p['s13'], p['dCP'], p['s14'], p['d14'], p['s15'], p['d15'],
            p['s24'], p['d24'], p['s25'], p['s34'], p['s35'], p['d35'],
            p['D21'], p['D31'], p['D41'], p['D51'], nubar=nubar)
    raise ValueError(d)


# ----------------------------------------------------------------------
# Density profiles.  Every one returns an ELECTRON NUMBER DENSITY [eV^3],
# so vcc_of() below must be called with is_number_density=True.
# ----------------------------------------------------------------------

L_SCALE = gd.L_SCALE_SUN
R_SUN = gd.SUN_RADIUS*gd.UNIT_KM
NE0 = gd.NUM_DENSITY_E_SUN_CENTRAL


def scalarize(y):
    """Unwrap a 0-d array to a plain scalar.

    The wrappers validate rho_func by calling it at L0 and requiring a float/int back; a
    profile written array-first returns a 0-d ndarray there and is rejected.  Every profile
    below returns `scalarize(...)` so it stays array-capable (keeping osc_prob's vectorized
    fast path) while still answering a single position with a number."""
    a = np.asarray(y)
    return a[()] if a.ndim == 0 else a


def solar_ne():
    """The package's own solar profile: Ne(r) = Ne(0) exp(-r/l_scale)."""
    return matter.exp_density_profile(NE0, L_SCALE)


def modulated_ne(amp=0.9, n_cycles=6.0, span=None):
    """Solar exponential modulated by a sine: crosses the resonance density repeatedly.
    This is the multi-resonance family the branch was developed against."""
    span = span if span is not None else 7.0*L_SCALE

    def ne(l):
        base = NE0*np.exp(-np.asarray(l, dtype=float)/L_SCALE)
        return scalarize(base*(1.0 + amp*np.sin(
            2.0*np.pi*n_cycles*np.asarray(l, dtype=float)/span)))
    return ne


def narrow_bump_ne(l_center, width, height_ratio=8.0, base_ratio=1.0e-3):
    """A quiet exponential background with ONE narrow Gaussian bump of the given width.
    Aimed squarely at the detector's fixed n_probe = 200 linear grid."""
    def ne(l):
        x = np.asarray(l, dtype=float)
        base = NE0*base_ratio*np.exp(-x/(4.0*L_SCALE))
        bump = NE0*base_ratio*height_ratio*np.exp(-0.5*((x - l_center)/width)**2)
        return scalarize(base + bump)
    return ne


def multi_bump_ne(centers, width, height_ratio=8.0, base_ratio=1.0e-3):
    centers = np.asarray(centers, dtype=float)

    def ne(l):
        x = np.asarray(l, dtype=float)[..., None]
        base = NE0*base_ratio*np.exp(-np.asarray(l, dtype=float)/(4.0*L_SCALE))
        bumps = NE0*base_ratio*height_ratio*np.exp(-0.5*((x - centers)/width)**2)
        return scalarize(base + bumps.sum(axis=-1))
    return ne


def sine_ne(period, amp=0.9, base_ratio=1.0e-2):
    """Pure sinusoidal modulation at a chosen spatial period -- for aliasing against n_probe."""
    def ne(l):
        x = np.asarray(l, dtype=float)
        return scalarize(NE0*base_ratio*(1.0 + amp*np.sin(2.0*np.pi*x/period)))
    return ne


def fourier_ne(rng, n_modes=8, span=None, base_ratio=1.0e-2, amp=0.8):
    """Random smooth profile with controlled bandwidth: a random Fourier sum, kept positive."""
    span = span if span is not None else 4.0*L_SCALE
    ks = rng.integers(1, n_modes + 1, size=n_modes)
    a = rng.normal(size=n_modes)
    b = rng.normal(size=n_modes)
    nrm = np.sqrt(np.sum(a**2 + b**2)) or 1.0
    a, b = a/nrm, b/nrm

    def ne(l):
        x = np.asarray(l, dtype=float)[..., None]
        ph = 2.0*np.pi*ks*x/span
        s = (a*np.cos(ph) + b*np.sin(ph)).sum(axis=-1)
        return scalarize(NE0*base_ratio*(1.0 + amp*s)
                         * np.exp(-np.asarray(l, dtype=float)/(6.0*L_SCALE)))
    return ne


def vcc_of(ne_func, nubar=False):
    """V_CC from an ELECTRON NUMBER DENSITY profile.

    TRAP GUARD: the 7th positional argument of vcc_func_from_rho_func is
    `density_is_of_number_of_electrons`, NOT `nubar`.  Passing nubar there silently
    yields a potential ~1e9 too small -- a vacuum reference that looks converged.
    Both are passed by KEYWORD here so the mistake is unrepresentable.
    """
    return matter.vcc_func_from_rho_func(
        ne_func, 0.0, 1.0, 0.5,
        nubar=nubar,
        density_matter_is_in_g_per_cm3=False,
        density_is_of_number_of_electrons=True)


# ----------------------------------------------------------------------
# Hamiltonian factory -- ONE parameter, always.
# ----------------------------------------------------------------------

def H_factory(d, osc_params, vcc_func, energy, nubar=False, h_vac=None):
    """Returns H_of_l(l), a genuine ONE-argument closure, built exactly as
    _osc_prob_hybrid_dispatch builds it: H(l) = (1/E) h_vac + VCC(l) * h_matt_proj."""
    h_vac = h_vac_for(d, osc_params, nubar=nubar) if h_vac is None else h_vac
    h_vac = np.asarray(h_vac, dtype=complex)
    proj = np.zeros((d, d), dtype=complex)
    proj[0, 0] = 1.0
    inv_E = 1.0/float(energy)

    def H_of_l(l):
        vcc = np.asarray(vcc_func(l))
        return inv_E*h_vac + vcc[..., None, None]*proj

    return H_of_l


# ----------------------------------------------------------------------
# Oracles
# ----------------------------------------------------------------------

def exact_U(H_func, l0, l1, dim, rtol=1e-12, atol=1e-14):
    """THE accuracy oracle: solve_ivp/DOP853, dU/dl = -i H(l) U.

    t_eval is left unset but we only ever read the endpoint; scipy still stores every
    accepted step, so keep l1-l0 modest at very low energy (documented trap)."""
    def rhs(l, y):
        return (-1j*np.asarray(H_func(l)) @ y.reshape(dim, dim)).ravel()

    sol = solve_ivp(rhs, (float(l0), float(l1)), np.eye(dim, dtype=complex).ravel(),
                    rtol=rtol, atol=atol, method='DOP853', dense_output=False)
    if not sol.success:
        raise RuntimeError('solve_ivp failed: ' + str(sol.message))
    return sol.y[:, -1].reshape(dim, dim)


def exact_U_many(H_func, l0, Ls, dim, rtol=1e-12, atol=1e-14):
    """Oracle at many baselines in one solve, via t_eval (which also suppresses the
    store-every-step behaviour that is ruinous at low energy)."""
    Ls = np.atleast_1d(np.asarray(Ls, dtype=float))

    def rhs(l, y):
        return (-1j*np.asarray(H_func(l)) @ y.reshape(dim, dim)).ravel()

    sol = solve_ivp(rhs, (float(l0), float(Ls[-1])), np.eye(dim, dtype=complex).ravel(),
                    rtol=rtol, atol=atol, method='DOP853', t_eval=Ls)
    if not sol.success:
        raise RuntimeError('solve_ivp failed: ' + str(sol.message))
    return np.array([sol.y[:, i].reshape(dim, dim) for i in range(len(Ls))])


def oracle_converged(H_func, l0, l1, dim, quoted_error):
    """Verify the oracle itself: tighten rtol to 1e-13 and confirm it moves by << the
    error being quoted.  Returns (moved_by, ok)."""
    Ua = exact_U(H_func, l0, l1, dim, rtol=1e-12, atol=1e-14)
    Ub = exact_U(H_func, l0, l1, dim, rtol=1e-13, atol=1e-15)
    moved = float(np.max(np.abs(P_of(Ua) - P_of(Ub))))
    return moved, moved < 0.1*max(quoted_error, 1e-15)


def P_of(U):
    """Probability matrix in the package's convention: P[i][f] = |U[f][i]|^2."""
    U = np.asarray(U)
    return np.transpose(U.real**2 + U.imag**2)


def unitarity(U):
    U = np.asarray(U)
    return float(np.max(np.abs(U.conj().T @ U - np.eye(U.shape[-1]))))


def maxabs(x):
    return float(np.max(np.abs(np.asarray(x))))


# ----------------------------------------------------------------------
# Warning capture
# ----------------------------------------------------------------------

class Caught:
    """Records every warning raised in the block, by category name."""

    def __init__(self):
        self.records = []

    def __enter__(self):
        self._cm = warnings.catch_warnings(record=True)
        self._w = self._cm.__enter__()
        warnings.simplefilter('always')
        return self

    def __exit__(self, *exc):
        self.records = list(self._w)
        self._cm.__exit__(*exc)
        return False

    @property
    def names(self):
        return sorted({r.category.__name__ for r in self.records})

    def has(self, name):
        return any(r.category.__name__ == name for r in self.records)


def fmt(x):
    return ('%.3e' % x) if isinstance(x, float) else str(x)
