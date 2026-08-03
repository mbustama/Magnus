# -*- coding: utf-8 -*-
"""Does the Battery 2 defect reach a USER, through the package's public entry points?

hybrid_propagator's certified=True is internal.  What matters is whether
osc_prob_matter_std_potential(..., strategy='auto') -- the default -- returns the wrong
number SILENTLY.  Under 'auto', certified=False makes the dispatcher return NotImplemented
and fall back to the general path; only certified=True is actually returned to the caller.

Also checks strategy='hybrid' (explicit) and strategy='magnus' (the escape hatch), and
whether the cumulative route (N >= 25) rescues the same profile.
"""

import numpy as np

import harness as H
import harness6 as H6
import magnus.oscprob as oscprob
from battery2 import L0, L1, E, bump_profile, ne_res_for

span = L1 - L0
p2 = H.params_for(2)
ner2 = ne_res_for(2, p2, E)
rng = np.random.default_rng(7)
lc = L0 + (0.37 + 0.2*rng.random())*span
TOL = 1e-3

CASES = [
    ('sub-threshold, w=1e-2 span', bump_profile(ner2, lc, 1e-2*span)),
    ('sub-threshold, w=3e-2 span', bump_profile(ner2, lc, 3e-2*span)),
    ('detection miss, w=1e-5 span', bump_profile(ner2, lc, 1e-5*span)),
    ('detection miss, w=3e-5 span', bump_profile(ner2, lc, 3e-5*span)),
    ('sinusoid, period = span/7',
     lambda l: ner2*(1.0 + 0.9*np.sin(2.0*np.pi*np.asarray(l, float)/(span/7.0)))),
]

print('%-30s %-9s %11s %-9s %-34s' % ('case', 'strategy', 'err', 'verdict', 'warnings'))
print('-'*100)

for name, ne in CASES:
    vcc = H.vcc_of(ne)
    H_of_l = H6.H_family('std', 2, E, vcc, p2)
    Pref = H.P_of(H.exact_U(H_of_l, L0, L1, 2))

    for strat in ('auto', 'hybrid', 'magnus'):
        with H.Caught() as c:
            P = np.asarray(oscprob.osc_prob_matter_std_potential(
                2, ne, E, L1, p2, L0=L0, density_is_of_number_of_electrons=True,
                strategy=strat))
        err = H.maxabs(P - Pref)
        silent = err > TOL and not c.names
        print('%-30s %-9s %11.3e %-9s %-34s'
              % (name, strat, err, 'SILENT' if silent else ('warned' if err > TOL else 'ok'),
                 ','.join(c.names) or '-'), flush=True)

    # And the cumulative route: the same profile as a 60-point baseline scan (N >= 25),
    # which HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS sends to the cumulative scan instead.
    Ls = np.linspace(0.05*L1, L1, 60)
    with H.Caught() as c:
        Pscan = np.asarray(oscprob.osc_prob_matter_std_potential(
            2, ne, E, Ls, p2, L0=L0, density_is_of_number_of_electrons=True,
            strategy='auto')).reshape(60, 2, 2)
    Uref = H.exact_U_many(H_of_l, L0, Ls, 2)
    Pr = np.array([H.P_of(U) for U in Uref])
    err = H.maxabs(Pscan - Pr)
    print('%-30s %-9s %11.3e %-9s %-34s'
          % (name, 'auto N=60', err,
             'SILENT' if (err > TOL and not c.names) else ('warned' if err > TOL else 'ok'),
             ','.join(c.names) or '-'), flush=True)
    print()
