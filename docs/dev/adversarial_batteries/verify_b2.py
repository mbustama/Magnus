# -*- coding: utf-8 -*-
"""Verify the Battery 2.1 failures are code bugs, not instrument bugs.

The worry: an ADAPTIVE oracle (DOP853) starting far from a very narrow bump could step
straight over it, exactly as the detector does -- which would make a spurious "error".
Three independent checks per failing configuration:

  1. rtol 1e-12 -> 1e-13 (the standard convergence check).
  2. max_step capped at width/10, forcing the solver to sample inside the bump.
  3. PIECEWISE solve: split the domain at the bump and compose U = U2 @ U1, so the solver
     is restarted with the feature at the centre of a short interval it cannot miss.
  4. An independent Magnus reference: osc_prob on a fixed grid with t_breakpoints ON the
     bump, refined until it stops moving.  (Cross-family check, not a Magnus-vs-Magnus
     accuracy claim: it only has to agree with solve_ivp.)
"""

import numpy as np

import harness as H
import magnus.adiabatic as ad
import magnus.oscprob as oscprob
from battery2 import L0, L1, E, bump_profile, ne_res_for
from scipy.integrate import solve_ivp

span = L1 - L0
p2 = H.params_for(2)
ner2 = ne_res_for(2, p2, E)
rng = np.random.default_rng(7)
lc = L0 + (0.37 + 0.2*rng.random())*span


def U_capped(H_of_l, l0, l1, dim, max_step):
    def rhs(l, y):
        return (-1j*np.asarray(H_of_l(l)) @ y.reshape(dim, dim)).ravel()
    sol = solve_ivp(rhs, (l0, l1), np.eye(dim, dtype=complex).ravel(),
                    rtol=1e-12, atol=1e-14, method='DOP853', max_step=max_step)
    assert sol.success
    return sol.y[:, -1].reshape(dim, dim)


def U_piecewise(H_of_l, l0, l1, dim, cuts):
    """Compose across cuts, restarting the solver at each -- the bump cannot be skipped."""
    edges = np.concatenate([[l0], np.sort(np.asarray(cuts, float)), [l1]])
    U = np.eye(dim, dtype=complex)
    for a, b in zip(edges[:-1], edges[1:]):
        if b <= a:
            continue
        U = H.exact_U(H_of_l, a, b, dim) @ U
    return U


print('%-10s %12s %12s %12s %12s %12s %12s'
      % ('width/span', 'hybrid err', 'oracle d(rtol)', 'oracle d(cap)', 'oracle d(pw)',
         'magnus ref d', 'VERDICT'))

for wf in [3e-2, 1e-2, 3e-5, 1e-5]:
    w = wf*span
    ne = bump_profile(ner2, lc, w)
    vcc = H.vcc_of(ne)
    H_of_l = H.H_factory(2, p2, vcc, E)

    U_hyb, win, cert = ad.hybrid_propagator(H_of_l, L0, L1, rtol=1e-3, atol=1e-3)
    U_a = H.exact_U(H_of_l, L0, L1, 2)                      # the standard oracle
    U_b = H.exact_U(H_of_l, L0, L1, 2, rtol=1e-13, atol=1e-15)
    U_c = U_capped(H_of_l, L0, L1, 2, max_step=w/10.0)
    U_d = U_piecewise(H_of_l, L0, L1, 2,
                      [lc - 8*w, lc - 2*w, lc, lc + 2*w, lc + 8*w])

    # Independent Magnus reference: fixed dense grid, breakpoints straddling the bump.
    bps = np.array([lc - 8*w, lc - 3*w, lc, lc + 3*w, lc + 8*w])
    bps = bps[(bps > L0) & (bps < L1)]
    P_m = None
    for n in (4000, 8000, 16000):
        Pn = oscprob.osc_prob(H_of_l, L0, L1, n_slabs=n, n_tpts_per_slab=100,
                              magnus_exp_order=6, rtol=None, atol=None,
                              t_breakpoints=bps)
        if P_m is not None and H.maxabs(np.asarray(Pn) - P_m) < 1e-9:
            P_m = np.asarray(Pn)
            break
        P_m = np.asarray(Pn)

    e_hyb = H.maxabs(H.P_of(U_hyb) - H.P_of(U_a))
    d_rtol = H.maxabs(H.P_of(U_a) - H.P_of(U_b))
    d_cap = H.maxabs(H.P_of(U_a) - H.P_of(U_c))
    d_pw = H.maxabs(H.P_of(U_a) - H.P_of(U_d))
    d_mag = H.maxabs(H.P_of(U_a) - P_m)

    verdict = 'CODE BUG' if (max(d_rtol, d_cap, d_pw, d_mag) < 0.05*e_hyb) else 'ORACLE SUSPECT'
    print('%-10.0e %12.3e %12.3e %12.3e %12.3e %12.3e   %s  (windows=%d cert=%s)'
          % (wf, e_hyb, d_rtol, d_cap, d_pw, d_mag, verdict, len(win), cert))
