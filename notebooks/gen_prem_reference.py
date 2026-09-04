# -*- coding: utf-8 -*-
r"""gen_prem_reference.py

An extended-precision reference for the PREM chord at cos(theta_z) = -0.9, the Earth
analogue of the exponential profile behind Figure 11.

WHY IT IS NOT THE SAME SCRIPT AS ``gen_mp_reference.py``.  Triple-Richardson removes the
N^-2, N^-4 and N^-6 terms of the midpoint rule, and those terms only exist if the
integrand is smooth across every slab.  The exponential profile is smooth everywhere;
PREM is not.  It is polynomial in radius *inside* each shell and discontinuous *between*
them, so a uniform ladder whose slabs straddle a boundary has no clean expansion for
Richardson to cancel -- and a ladder can agree with itself beautifully while being wrong,
which is not hypothetical.  On the smooth profile, at 4nu, bases 32-256 converged onto
0.1926988 with a self-consistency of 5.9e-16 and the true answer was 0.1926925: wrong by
6.3e-06, certified to machine precision.  Only an independent method caught it.

SO THE SLABS ARE SEGMENT-ALIGNED.  The chord's sixteen layer crossings are always slab
boundaries.  Refinement multiplies the slab count *within* each segment, so every slab
lies inside one shell, the midpoint rule keeps its even-power expansion there, and
Richardson cancels what it is built to cancel.  This is the same thing Magnus itself does
with ``t_breakpoints``, for the same reason.

TWO ACCEPTANCE CRITERIA, NOT ONE.  Self-convergence is necessary and demonstrably not
sufficient, so a point is accepted only when the ladder has converged *and* an
independent double-precision DOP853 integration agrees with it to ``ODE_GATE``.

    python notebooks/gen_prem_reference.py

Writes ``prem_chord_reference.json``, checkpointing after every energy.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import json
import pathlib
import platform
import sys
import time

import numpy as np
from mpmath import mp
from scipy.integrate import solve_ivp

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import gen_profile_benchmarks as gpb                      # noqa: E402
import prem_chord_common as pcc                           # noqa: E402

from magnus import globaldefs as gd                       # noqa: E402

HERE = pathlib.Path(__file__).parent
OUT = HERE/'prem_chord_reference.json'

DPS = 30
TARGET = 1.0e-15
# Every point of the first run hit the old cap of 256 without reaching TARGET, leaving
# self-convergence at 2.3e-11 (3nu) to 6.0e-09 (4nu) -- which would become the new floor
# the moment the methods reach past it.  Triple Richardson converges as N^-8, so closing
# a factor of 6e6 needs roughly six times the slabs; 2048 leaves margin.  MULT0 starts
# where the old run ended, since everything below it is known to be insufficient and the
# ladders are memoized anyway.
MULT0 = 128
MULT_CAP = 2048
ODE_GATE = 1.0e-8
ENERGIES = gpb.exponential_profile()['energies']          # the same twelve, for analogy


def segment_slabs(edges, mult):
    """Slab edges that always include the layer boundaries.

    Each segment gets a slab count proportional to its length, never fewer than one, so
    refinement is uniform in the slab *width* rather than in the count per segment --
    the chord's segments differ by three orders of magnitude in length.
    """
    seg = np.diff(edges)
    n = np.maximum(1, np.round(mult*seg/seg.max()).astype(int))
    out = [np.linspace(edges[i], edges[i + 1], n[i] + 1)[:-1] for i in range(len(seg))]
    return np.concatenate(out + [edges[-1:]])


def mp_product(hf, cuts, d):
    """Midpoint slab product over the given slab edges, carried in mpmath."""
    mids = 0.5*(cuts[:-1] + cuts[1:])
    widths = np.diff(cuts)
    Hs = hf(mids)
    U = mp.eye(d)
    for k in range(len(widths)):
        M = mp.matrix([[mp.mpc(complex(Hs[k][i, j])) for j in range(d)]
                       for i in range(d)])
        U = mp.expm(-1j*M*mp.mpf(float(widths[k])))*U
    return U


def richardson3(P, ms):
    Q = {m: (4*P[2*m] - P[m])/3 for m in ms[:-1]}
    R = {m: (16*Q[2*m] - Q[m])/15 for m in ms[:-2]}
    T = (64*R[ms[1]] - R[ms[0]])/63
    return T, float(max(abs(x) for x in (T - R[ms[1]])))


def ode_reference(hf, L, d, nu_i, nu_f, rtol=1.0e-12):
    """Independent DOP853 integration of dU/dl = -i H U, in double precision."""
    def rhs(x, y):
        U = y.reshape(d, d, 2)
        dU = -1j*np.asarray(hf(x)) @ (U[..., 0] + 1j*U[..., 1])
        return np.stack([dU.real, dU.imag], axis=-1).ravel()

    y0 = np.stack([np.eye(d), np.zeros((d, d))], axis=-1).ravel()
    sol = solve_ivp(rhs, (0.0, L), y0, method='DOP853', rtol=rtol, atol=rtol*1.0e-2)
    U = sol.y[:, -1].reshape(d, d, 2)
    return abs((U[..., 0] + 1j*U[..., 1])[nu_f, nu_i])**2


def reference_at(hf, edges, d, nu_i, nu_f, ode):
    P, mult = {}, MULT0
    while True:
        ms = (mult, 2*mult, 4*mult, 8*mult)
        for m in ms:
            if m not in P:
                U = mp_product(hf, segment_slabs(edges, m), d)
                P[m] = mp.matrix([[abs(U[i, j])**2 for j in range(d)]
                                  for i in range(d)])
        T, sc = richardson3(P, ms)
        val = float(T[nu_f, nu_i])
        if (sc < TARGET and abs(val - ode) < ODE_GATE) or mult >= MULT_CAP:
            return val, sc, mult, len(segment_slabs(edges, 8*mult)) - 1
        mult *= 2


def main():
    ch = pcc.chord()
    L, edges = ch['baseline'], ch['edges']
    mp.dps = DPS
    out = {'what': 'segment-aligned triple-Richardson mpmath reference, PREM chord',
           'costhz': ch['costhz'], 'baseline_km': ch['baseline_km'],
           'baseline_inv_ev': L, 'layer_edges_inv_ev': [float(x) for x in edges],
           'energy_ev': [float(e) for e in ENERGIES], 'dps': DPS,
           'target_self_convergence': TARGET, 'machine': platform.platform(),
           'cases': []}
    if OUT.exists():
        out = json.loads(OUT.read_text())

    for d in (2, 3, 4, 5):
        case = next((c for c in out['cases'] if c['flavours'] == d), None)
        if case is None:
            case = {'flavours': d, 'reference': [], 'self_convergence': [],
                    'mult': [], 'n_slabs': [], 'ode_crosscheck': []}
            out['cases'].append(case)
        print('--- %d flavours ---' % d, file=sys.stderr, flush=True)
        for j, e in enumerate(ENERGIES):
            if j < len(case['reference']):
                continue
            t0 = time.time()
            hf = gpb.h_of(d, ch['vcc'], e)
            ode = ode_reference(hf, L, d, gd.NUMU, gd.NUMU)
            val, sc, mult, ns = reference_at(hf, edges, d, gd.NUMU, gd.NUMU, ode)
            case['reference'].append(val)
            case['self_convergence'].append(sc)
            case['mult'].append(mult)
            case['n_slabs'].append(ns)
            case['ode_crosscheck'].append(abs(val - ode))
            OUT.write_text(json.dumps(out, indent=1))
            print('  E %2d/%d  %5.2f GeV  P=%.15f  self-conv %.1e  mult %4d '
                  '(%d slabs)  vs DOP853 %.1e  (%.0fs)'
                  % (j + 1, len(ENERGIES), e/gd.UNIT_GEV, val, sc, mult, ns,
                     abs(val - ode), time.time() - t0), file=sys.stderr, flush=True)
    print('wrote %s' % OUT.name, file=sys.stderr)


if __name__ == '__main__':
    main()
