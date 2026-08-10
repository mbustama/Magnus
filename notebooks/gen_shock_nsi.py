# -*- coding: utf-8 -*-
r"""gen_shock_nsi.py

The supernova shock of notebook 14 with non-standard interactions, for Mag(nu)s and
NuOscProbExact.  Writes ``external_shock_nsi.json``, which notebook 25 reads.

    python notebooks/gen_shock_nsi.py > notebooks/external_shock_nsi.json

THE CONVENTIONS WERE CHECKED FIRST, AND THEY MATCH EXACTLY.  Both codes take six
dimensionless ``eps`` and both build the matter term as ``V_CC * (diag(1,0,0) + eps)`` --
NuOscProbExact writes the standard piece into the matrix as ``1 + eps_ee`` while Mag(nu)s
adds it separately, which looks like an off-by-one in the ee entry and is not.  Measured on
constant density (3 g/cm^3, 1300 km, 2 GeV, ``eps_ee = 0.10``, ``eps_em = 0.05``), the same
``eps`` handed to both gives P(numu -> nue) agreeing to **1.7e-16**, and shifting either
convention by one puts them 3.7e-02 apart.  That is the check this repository's own rule
demands before any cross-code ratio is quoted, and unlike the V_CC comparison in section 6
it comes out clean: there is no offset to correct here.

REFEREE: an adaptive DOP853 integration of the same Hamiltonian, which is neither a Magnus
expansion nor a slab product.  ``shock_reference.json`` cannot serve -- it freezes the
*standard* three-flavour problem -- so a reference is computed here.  The accumulated phase
is set by D31 as in the standard case, about 1.5e4 radians, so the referee costs what
notebook 14's did rather than what the 3+1 case would (see ``gen_shock_4nu.py`` on why an
eV-scale splitting cannot be refereed at all).

DRIVERS.  Mag(nu)s in one call for all 61 baselines with the fronts declared through
``t_breakpoints``; the cumulative scan engages by itself at this many baselines, verified by
``cumulative=True`` changing the answer by exactly zero.  NuOscProbExact by composing
``evolution_operator_3nu_slabs`` across legs cut at the declared fronts and at every target,
with a per-front slab floor so the fronts are resolved rather than jumped -- the allocation
that, left proportional to length, made the code look like it floored at 5e-07.
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
import magnus.matter as matter                                    # noqa: E402
import magnus.oscprob as oscprob                                  # noqa: E402

NOTEBOOK = '14_magnus_supernova_shock.ipynb'
WIDTH = 1.0e-3
# Large enough to be legible against the standard curve and still a perturbation.  eps_em
# is the one that moves nu_e <-> nu_mu directly, which is the channel plotted.
EPS = dict(eps_ee=0.15, eps_em=0.05, eps_et=0.0, eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)
MG_DIALS = (2000, 8000, 32000, 128000)
NPE_DIALS = (8192, 32768, 131072, 524288)
MAX_STACK_BYTES = 192*1024**2


def have(name):
    try:
        __import__(name)
        return True
    except Exception:
        return False


HAVE_NPE = have('slabs')


def notebook_namespace():
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
H_VAC3 = np.asarray(NS['hvac3'])

EPS_MATRIX = np.array([[1.0 + EPS['eps_ee'], EPS['eps_em'], EPS['eps_et']],
                       [np.conj(EPS['eps_em']), EPS['eps_mm'], EPS['eps_mt']],
                       [np.conj(EPS['eps_et']), np.conj(EPS['eps_mt']), EPS['eps_tt']]],
                      dtype=complex)


def vcc_of(width_frac):
    """V_CC along the ray, from notebook 14's own builder, checked independently."""
    h_fixed = NS['make_H'](NS['sn_shock_ne'](width_frac))
    hvac_over_e = H_VAC3/ENERGY

    def f(l):
        return np.asarray(h_fixed(l))[..., 0, 0] - hvac_over_e[0, 0]

    check = matter.vcc_func_from_rho_func(
        NS['sn_shock_ne'](width_frac), 0.0, 1.0, 0.5, nubar=False,
        density_matter_is_in_g_per_cm3=False, density_is_of_number_of_electrons=True)
    probe = np.linspace(L0, L1, 97)
    got, want = np.asarray(f(probe)), np.asarray(check(probe))
    if not np.allclose(got, want, rtol=1.0e-12, atol=0.0):
        raise SystemExit('V_CC extracted from notebook 14 disagrees with an independent '
                         'construction by up to %.2e' % float(np.max(np.abs(got - want))))
    return f


def h_nsi_of(width_frac, energy):
    """H = h_vac/E + V_CC (diag(1,0,0) + eps).  One definition, given to every code."""
    vcc = vcc_of(width_frac)

    def f(l):
        la = np.atleast_1d(np.asarray(l, dtype=float))
        out = np.broadcast_to((H_VAC3/energy).astype(complex), (len(la), 3, 3)).copy()
        out += np.asarray(vcc(la))[:, None, None]*EPS_MATRIX
        return out[0] if np.ndim(l) == 0 else out
    return f


def referee(width_frac):
    h_of = h_nsi_of(width_frac, ENERGY)

    def rhs(l, y):
        return (-1j*np.asarray(h_of(l)) @ y.reshape(3, 3)).ravel()

    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (L0, float(Ls[-1])), np.eye(3, dtype=complex).ravel(),
                    rtol=1.0e-12, atol=1.0e-14, method='DOP853', t_eval=Ls)
    if not sol.success:
        raise SystemExit('DOP853 failed: %s' % sol.message)
    u = np.array([sol.y[:, i].reshape(3, 3) for i in range(len(Ls))])
    p = np.swapaxes(u.real**2 + u.imag**2, -1, -2)
    print('  referee: %.1f s, unitarity %.2e'
          % (time.perf_counter() - t0, float(np.max(np.abs(p.sum(axis=2) - 1.0)))),
          file=sys.stderr)
    return p


def timed(call, repeat=3, min_block=0.05, budget=8.0):
    t0 = time.perf_counter()
    call()
    if time.perf_counter() - t0 > budget:
        return time.perf_counter() - t0
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


def magnus_along_ray(width_frac, n_slabs, nsi=True):
    fn = oscprob.osc_prob_matter_nsi if nsi else oscprob.osc_prob_matter_std_potential
    args = ((3, NS['sn_shock_ne'](width_frac), ENERGY, Ls, PARAMS3, EPS) if nsi else
            (3, NS['sn_shock_ne'](width_frac), ENERGY, Ls, PARAMS3))
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.asarray(fn(
            *args, L0=L0, density_is_of_number_of_electrons=True,
            t_breakpoints=NS['shock_breakpoints'](width_frac),
            n_slabs=n_slabs, max_n_slabs=4*n_slabs,
            rtol=1.0e-12, atol=1.0e-14)).reshape(len(Ls), 3, 3)


def npe_along_ray(width_frac, n_slabs):
    import slabs as npe

    h_of = h_nsi_of(width_frac, ENERGY)
    bps = np.asarray(NS['shock_breakpoints'](width_frac), dtype=float)
    h = (L1 - L0)/float(n_slabs)
    floor = max(8, n_slabs//1000)
    edges = np.unique(np.concatenate([bps[bps < Ls[0]], [L0, Ls[0]]]))
    edges = edges[(edges >= L0) & (edges <= Ls[0])]

    used = 0
    u = np.eye(3, dtype=complex)
    step = max(1, int(MAX_STACK_BYTES//(9*16)))
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = max(floor, int(round((hi - lo)/h)))
        e = np.linspace(lo, hi, m + 1)
        for j0 in range(0, m, step):
            sub = e[j0:j0 + step + 1]
            u = np.asarray(npe.evolution_operator_3nu_slabs(
                np.asarray(h_of(0.5*(sub[:-1] + sub[1:])), dtype=complex),
                np.diff(sub))) @ u
        used += m

    d_l = Ls[1] - Ls[0]
    m = max(1, int(round(d_l/h)))
    e = np.linspace(0.0, d_l, m + 1)
    mid = Ls[:-1, None] + 0.5*(e[:-1] + e[1:])[None, :]
    stack = np.asarray(h_of(mid.ravel()), dtype=complex).reshape(len(Ls) - 1, m, 3, 3)
    u_batch = np.asarray(npe.evolution_operator_3nu_slabs(stack, np.diff(e)))
    used += (len(Ls) - 1)*m

    ops = [u]
    for j in range(len(Ls) - 1):
        u = u_batch[j] @ u
        ops.append(u)
    ops = np.array(ops)
    return np.swapaxes(ops.real**2 + ops.imag**2, -1, -2), used


def main():
    out = {'note': ('supernova shock of notebook 14 with NSI; produced by '
                    'notebooks/gen_shock_nsi.py'),
           'machine': platform.platform(),
           'width': WIDTH, 'energy_mev': ENERGY/gd.UNIT_MEV,
           'eps': {k: float(v) for k, v in EPS.items()},
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
    out['magnus_P_ee'] = [float(x)
                          for x in magnus_along_ray(WIDTH, MG_DIALS[-1])[:, 0, 0]]
    # The standard three-flavour curve, so the NSI departure has something to be read
    # against.  A BSM curve on its own says nothing about size.
    out['standard_P_ee'] = [float(x) for x in
                            magnus_along_ray(WIDTH, MG_DIALS[-1], nsi=False)[:, 0, 0]]

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
