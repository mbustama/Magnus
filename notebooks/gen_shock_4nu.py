# -*- coding: utf-8 -*-
r"""gen_shock_4nu.py

The supernova shock of notebook 14, at 3+1, for Mag(nu)s and NuOscProbExact.  Writes
``external_shock_4nu.json``, which notebook 25 reads.

    python notebooks/gen_shock_4nu.py > notebooks/external_shock_4nu.json

WHY THIS IS A SEPARATE CASE AND NOT A FLAG ON THE 3-FLAVOUR ONE.  ``shock_reference.json``
freezes a DOP853 solution of the *three-flavour* Schroedinger equation, so it cannot referee
four flavours; a reference is computed here instead.  Everything else is shared: the profile,
the ray, the energy and the sampled baselines all come from notebook 14's own cells, so there
is one definition of the shock and both notebooks read it.

THE STERILE STATE FEELS THE MEDIUM, and getting that wrong is invisible.  The matter term is
``diag(1, 0, 0, r/2) * V_CC`` with ``r = n_n/n_p``: actives share V_NC and it cancels, a sterile
state feels neither current and keeps ``-V_NC = (r/2) V_CC``.  Omitting it costs 0.29 in
probability on a PREM chord, flat in tolerance, so no refinement reveals it -- and notebook 25's
own PREM referee carried the omission until it was found by a referee disagreeing with every
code it refereed.  The projector here therefore comes from
``matter.matter_potential_projector`` and is never written out.

DRIVERS, as in ``gen_shock_benchmarks.py`` and for the same reasons: Mag(nu)s cumulatively in
one call for all 61 baselines, NuOscProbExact by composing ``evolution_operator_4nu_slabs``
across legs cut at the declared fronts and at every target, with a per-front slab floor so the
0.07 km fronts are resolved.  Both are told where the front is.  Both are chunked to a memory
budget, because the unguarded phrasing of this asks for gigabytes.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import contextlib
import json
import pathlib
import platform
import sys
import time
import warnings

import numpy as np
from scipy.integrate import solve_ivp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'src'))
sys.path.insert(0, str(HERE))

import magnus.globaldefs as gd                                    # noqa: E402
import magnus.hamiltonians as hamiltonians                        # noqa: E402
import magnus.matter as matter                                    # noqa: E402
import magnus.oscprob as oscprob                                  # noqa: E402

NOTEBOOK = '14_magnus_supernova_shock.ipynb'
WIDTH = 1.0e-3                    # the simulation-smeared front; see the note in main()
# D41 IS 1e-2 AND NOT THE eV-SCALE VALUE THE PREM SECTION USES, and the reason is the
# referee rather than the physics.  Accumulated phase over this 70 000 km ray at 15 MeV:
#
#     D31 = 2.5e-3    1.5e+04 rad       2 352 oscillation lengths   (the 3nu case)
#     D41 = 1.0       5.9e+06 rad     940 981 oscillation lengths
#     D41 = 1e-2      5.9e+04 rad       9 410 oscillation lengths
#
# An adaptive DOP853 reference has to resolve every one of those oscillations, so at an
# eV-scale splitting it needs of order a day for a single width -- measured by starting
# one and killing it after an hour with the referee still unfinished.  At 1e-2 the
# referee costs about four times the three-flavour one, which is affordable, and the
# sterile oscillation is still fast against the shock structure, which is what the case
# is about.  An eV-scale sterile state on this ray is not a case that can be *refereed*
# at all by an independent integrator, and that is worth knowing before trying.
STERILE = dict(s14=np.sqrt(0.10), s24=np.sqrt(0.10), s34=0.0, d14=0.0, d24=0.0,
               D41=1.0e-2)
MG_DIALS = (2000, 8000, 32000)
NPE_DIALS = (8192, 32768, 131072)
MAX_STACK_BYTES = 192*1024**2
BYTES_PER_SLAB_PER_ENERGY = 16*16                                 # 4x4 complex128


def have(name):
    try:
        __import__(name)
        return True
    except Exception:
        return False


HAVE_NPE = have('slabs')


def notebook_namespace():
    """Notebook 14's own definition cells; its stdout is kept off this script's."""
    import matplotlib
    matplotlib.use('Agg')
    import make_notebooks

    ns = {'__name__': '__notebook__'}
    with contextlib.redirect_stdout(sys.stderr):
        for cell in make_notebooks.books[NOTEBOOK].cells:
            if cell.cell_type != 'code':
                continue
            if 'def measure(' in cell.source:
                break
            exec(compile(cell.source, '<%s>' % NOTEBOOK, 'exec'), ns)
    return ns


NS = notebook_namespace()
L0, L1 = NS['L0'], NS['L1']
Ls = np.asarray(NS['Ls'], dtype=float)
ENERGY, PARAMS3 = NS['ENERGY'], NS['params3']
OSC4 = dict(PARAMS3, **STERILE)

H_VAC4 = np.asarray(hamiltonians.hamiltonian_4nu_vacuum_energy_independent(
    PARAMS3['s12'], PARAMS3['s23'], PARAMS3['s13'], PARAMS3['dCP'],
    STERILE['s14'], 0.0, STERILE['s24'], 0.0, STERILE['s34'],
    PARAMS3['D21'], PARAMS3['D31'], STERILE['D41']))
PROJ4 = matter.matter_potential_projector(4)


def vcc_of(width_frac):
    r"""V_CC along the ray, taken from notebook 14's own builder.

    ``make_H`` returns ``h_vac3/ENERGY + V_CC(l) * diag(1,0,0)``, so subtracting the
    constant vacuum part leaves ``V_CC`` in the (0,0) entry exactly.  It is asserted
    against an independently constructed potential rather than trusted, because a
    silently wrong potential is precisely the failure this whole exercise keeps finding.
    """
    h_fixed = NS['make_H'](NS['sn_shock_ne'](width_frac))
    hvac3 = np.asarray(NS['hvac3'])/ENERGY

    def f(l):
        return np.asarray(h_fixed(l))[..., 0, 0] - hvac3[0, 0]

    check = matter.vcc_func_from_rho_func(
        NS['sn_shock_ne'](width_frac), 0.0, 1.0, 0.5, nubar=False,
        density_matter_is_in_g_per_cm3=False, density_is_of_number_of_electrons=True)
    probe = np.linspace(L0, L1, 97)
    got, want = np.asarray(f(probe)), np.asarray(check(probe))
    if not np.allclose(got, want, rtol=1.0e-12, atol=0.0):
        raise SystemExit('V_CC extracted from notebook 14 disagrees with an independent '
                         'construction by up to %.2e -- the notebook has changed'
                         % float(np.max(np.abs(got - want))))
    return f


def h4_of(width_frac, energy):
    """The one 4x4 Hamiltonian every code in this file is given."""
    vcc = vcc_of(width_frac)

    def f(l):
        la = np.atleast_1d(np.asarray(l, dtype=float))
        out = np.broadcast_to(H_VAC4/energy, (len(la), 4, 4)).copy()
        out += np.asarray(vcc(la))[:, None, None]*PROJ4
        return out[0] if np.ndim(l) == 0 else out
    return f


def referee(width_frac):
    """DOP853 over the ray, evaluated at the 61 targets.  Neither code's method."""
    h_of = h4_of(width_frac, ENERGY)

    def rhs(l, y):
        return (-1j*np.asarray(h_of(l)) @ y.reshape(4, 4)).ravel()

    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (L0, float(Ls[-1])), np.eye(4, dtype=complex).ravel(),
                    rtol=1.0e-12, atol=1.0e-14, method='DOP853', t_eval=Ls)
    if not sol.success:
        raise SystemExit('DOP853 failed: %s' % sol.message)
    u = np.array([sol.y[:, i].reshape(4, 4) for i in range(len(Ls))])
    p = np.swapaxes(u.real**2 + u.imag**2, -1, -2)
    print('  referee: %.1f s, unitarity %.2e'
          % (time.perf_counter() - t0, float(np.max(np.abs(p.sum(axis=2) - 1.0)))),
          file=sys.stderr)
    return p


def timed(call, repeat=3, min_block=0.05, budget=8.0):
    t0 = time.perf_counter()
    call()
    first = time.perf_counter() - t0
    if first > budget:
        return first
    reps = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        el = time.perf_counter() - t0
        if el >= min_block:
            break
        reps *= 2
    best = el/reps
    for _ in range(repeat - 1):
        if best*reps > budget:
            break
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        best = min(best, (time.perf_counter() - t0)/reps)
    return best


def magnus_along_ray(width_frac, n_slabs):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.asarray(oscprob.osc_prob_matter_std_potential(
            4, NS['sn_shock_ne'](width_frac), ENERGY, Ls, OSC4, L0=L0,
            density_is_of_number_of_electrons=True,
            t_breakpoints=NS['shock_breakpoints'](width_frac),
            cumulative=True, n_slabs=n_slabs, max_n_slabs=4*n_slabs,
            rtol=1.0e-12, atol=1.0e-14)).reshape(len(Ls), 4, 4)


def npe_along_ray(width_frac, n_slabs):
    """The 61 probabilities, by composing NuOscProbExact 4nu slab-chain operators."""
    import slabs as npe

    h_of = h4_of(width_frac, ENERGY)
    bps = np.asarray(NS['shock_breakpoints'](width_frac), dtype=float)
    h = (L1 - L0)/float(n_slabs)
    floor = max(8, n_slabs//1000)
    edges = np.unique(np.concatenate([bps[bps < Ls[0]], [L0, Ls[0]]]))
    edges = edges[(edges >= L0) & (edges <= Ls[0])]

    used = 0
    u = np.eye(4, dtype=complex)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = max(floor, int(round((hi - lo)/h)))
        e = np.linspace(lo, hi, m + 1)
        step = max(1, int(MAX_STACK_BYTES//BYTES_PER_SLAB_PER_ENERGY))
        for j0 in range(0, m, step):
            sub = e[j0:j0 + step + 1]
            u = np.asarray(npe.evolution_operator_4nu_slabs(
                np.asarray(h_of(0.5*(sub[:-1] + sub[1:])), dtype=complex),
                np.diff(sub))) @ u
        used += m

    d_l = Ls[1] - Ls[0]
    m = max(1, int(round(d_l/h)))
    e = np.linspace(0.0, d_l, m + 1)
    mid = Ls[:-1, None] + 0.5*(e[:-1] + e[1:])[None, :]
    stack = np.asarray(h_of(mid.ravel()), dtype=complex).reshape(len(Ls) - 1, m, 4, 4)
    u_batch = np.asarray(npe.evolution_operator_4nu_slabs(stack, np.diff(e)))
    used += (len(Ls) - 1)*m

    ops = [u]
    for j in range(len(Ls) - 1):
        u = u_batch[j] @ u
        ops.append(u)
    ops = np.array(ops)
    return np.swapaxes(ops.real**2 + ops.imag**2, -1, -2), used


def main():
    out = {'note': ('supernova shock of notebook 14 at 3+1; produced by '
                    'notebooks/gen_shock_4nu.py'),
           'machine': platform.platform(),
           'width': WIDTH, 'energy_mev': ENERGY/gd.UNIT_MEV,
           'sterile': {k: float(v) for k, v in STERILE.items()},
           'n_targets': len(Ls),
           'targets_km': [float(x)/gd.UNIT_KM for x in Ls],
           'series': []}

    ref = referee(WIDTH)
    out['reference_unitarity'] = float(np.max(np.abs(ref.sum(axis=2) - 1.0)))
    out['reference_P_ee'] = [float(x) for x in ref[:, 0, 0]]

    mg = []
    for n in MG_DIALS:
        p = magnus_along_ray(WIDTH, n)
        t = timed(lambda n=n: magnus_along_ray(WIDTH, n))
        mg.append(dict(label=str(n), n_slabs=n,
                       us_per_probability=1.0e6*t/len(Ls),
                       max_abs_error=float(np.max(np.abs(p - ref)))))
        print('  magnus n=%-8d %.3e  %.0f us/prob'
              % (n, mg[-1]['max_abs_error'], mg[-1]['us_per_probability']),
              file=sys.stderr)
    out['series'].append({'name': 'Magnus', 'dial': 'n_slabs', 'points': mg})
    out['magnus_P_ee'] = [float(x) for x in magnus_along_ray(WIDTH, MG_DIALS[-1])[:, 0, 0]]

    if HAVE_NPE:
        pts = []
        for n in NPE_DIALS:
            p, used = npe_along_ray(WIDTH, n)
            t = timed(lambda n=n: npe_along_ray(WIDTH, n))
            pts.append(dict(label=str(n), n_slabs=n, n_slabs_used=used,
                            us_per_probability=1.0e6*t/len(Ls),
                            max_abs_error=float(np.max(np.abs(p - ref)))))
            print('  npe    n=%-8d %.3e  %.0f us/prob'
                  % (n, pts[-1]['max_abs_error'], pts[-1]['us_per_probability']),
                  file=sys.stderr)
        out['series'].append({'name': 'NuOscProbExact', 'dial': 'n_slabs',
                              'points': pts})
        out['npe_P_ee'] = [float(x)
                           for x in npe_along_ray(WIDTH, NPE_DIALS[-1])[0][:, 0, 0]]

    json.dump(out, sys.stdout, indent=1)


if __name__ == '__main__':
    main()
