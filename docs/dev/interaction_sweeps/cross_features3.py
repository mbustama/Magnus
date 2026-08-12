"""Third interaction sweep: cumulative, parallel workers, breakpoints, layer overrides.

Four routes that must not change the answer, and one guard that must not have a hole:

  * `cumulative` -- one propagation returning every baseline, against one call per baseline.
    Notebook 25 records that driving Magnus point-by-point instead cost a wrong answer, so
    the two routes existing means they have to agree.
  * `n_jobs` -- loky re-imports the package inside each worker, and a module-level global
    does not survive a process boundary.  There is already a test for EXPM_BACKEND not
    surviving it; anything else module-level is suspect for the same reason.
  * `t_breakpoints` -- the Earth wrappers place them on the PREM shell crossings for you.
    Passing the same set by hand has to be a no-op.
  * the per-layer `electron_fraction_*` overrides against the sterile projector guard added
    in O50, which computes its target from whichever composition is actually in force.
"""
import sys
import warnings

import numpy as np

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.oscprob as op

P = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
STD = {'s12': P['s12'], 's23': P['s23'], 's13': P['s13'], 'dCP': P['dCP']}
D3 = {'D21': P['D21'], 'D31': P['D31']}
COSTHZ = -0.7
CHORD = earth.distance_traveled_inside_earth(COSTHZ)
L_EARTH = CHORD*gd.CONV_KM_TO_INV_EV
E = 1.0*gd.UNIT_GEV


def worst(a, b):
    return float(np.max(np.abs(np.asarray(a, dtype=float) - np.asarray(b, dtype=float))))


def main():
    bad = []

    print('=== cumulative: one propagation for many baselines vs one call each ===')
    Ls = np.linspace(0.2*L_EARTH, L_EARTH, 9)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        cum = np.asarray(op.osc_prob_3nu_matter_constant_density(
            E, Ls, 2.8, **STD, **D3, density_matter_is_in_g_per_cm3=True,
            cumulative=True, rtol=1e-8, atol=1e-10), dtype=float)
        one = np.array([np.asarray(op.osc_prob_3nu_matter_constant_density(
            E, float(L), 2.8, **STD, **D3, density_matter_is_in_g_per_cm3=True,
            rtol=1e-8, atol=1e-10), dtype=float) for L in Ls])
    w = worst(cum, one)
    ok = w < 1e-6
    bad += [] if ok else [('cumulative', w)]
    print('  constant density, 9 baselines: worst |dP| = %.3e %s'
          % (w, '' if ok else '  <-- DISAGREE'))

    print()
    print('=== n_jobs: parallel workers against serial ===')
    energies = np.logspace(0.0, 1.0, 8)*gd.UNIT_GEV
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        # at a tolerance where the ladder has actually converged.  At the default
        # rtol=1e-3 the two differ by 1.2e-03 -- within what was asked for, but not
        # bitwise, because splitting the slabs changes the arithmetic order and the
        # stopping test compares successive levels.  That is documented, not a defect;
        # demanding bitwise agreement here was the harness being wrong about the contract.
        serial = op.osc_prob_3nu_earth(energies, costhz=COSTHZ, L=L_EARTH,
                                       **STD, **D3, n_jobs=1, rtol=1e-9, atol=1e-11)
        par = op.osc_prob_3nu_earth(energies, costhz=COSTHZ, L=L_EARTH,
                                    **STD, **D3, n_jobs=2, rtol=1e-9, atol=1e-11)
    w = worst(serial, par)
    ok = w < 1e-9
    bad += [] if ok else [('n_jobs', w)]
    print('  3nu earth, 8 energies: worst |dP| = %.3e %s'
          % (w, '' if ok else '  <-- DISAGREE'))

    # the same question for the newly added keyword: does `angles` survive a worker?
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        a = op.osc_prob_3nu_earth(energies, costhz=COSTHZ, L=L_EARTH,
                                  **STD, **D3, n_jobs=2)
        th = {k: float(np.degrees(np.arcsin(P[k]))) for k in ('s12', 's23', 's13')}
        b = op.osc_prob_3nu_earth(energies, costhz=COSTHZ, L=L_EARTH,
                                  s12=th['s12'], s23=th['s23'], s13=th['s13'],
                                  dCP=float(np.degrees(P['dCP'])), **D3,
                                  angles='deg', n_jobs=2)
    w = worst(a, b)
    ok = w < 1e-9
    bad += [] if ok else [('angles across workers', w)]
    print('  angles=deg across 2 workers: worst |dP| = %.3e %s'
          % (w, '' if ok else '  <-- DISAGREE'))

    print()
    print('=== t_breakpoints: passing the PREM crossings by hand is a no-op ===')
    edges = earth.prem_layer_edges_along_chord(COSTHZ)*gd.UNIT_KM
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        auto = op.osc_prob_3nu_earth(E, costhz=COSTHZ, L=L_EARTH, **STD, **D3)
        hand = op.osc_prob_3nu_earth(E, costhz=COSTHZ, L=L_EARTH, **STD, **D3,
                                     t_breakpoints=edges)
    w = worst(auto, hand)
    ok = w < 1e-12
    bad += [] if ok else [('t_breakpoints', w)]
    print('  3nu earth: worst |dP| = %.3e %s' % (w, '' if ok else '  <-- DISAGREE'))

    print()
    print('=== the O50 guard sees the per-layer overrides, not just the scalar ===')
    ster = dict(s14=0.15, d14=0.0, s24=0.10, d24=0.0, s34=0.05, D41=1.0)
    def warned(**kw):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            op.osc_prob_4nu_earth(E, costhz=-0.95,
                                  L=earth.distance_traveled_inside_earth(-0.95)*gd.CONV_KM_TO_INV_EV,
                                  **STD, **D3, **ster, **kw)
        return [c for c in caught
                if c.category is gd.SterileMatterCompositionWarning]
    cases = [
        ('layered defaults, r=1.0', {}, True),
        ('all layers pinned to Y_e=0.5, r=1.0', dict(
            electron_fraction_core=0.5, electron_fraction_mantle=0.5,
            electron_fraction_crust=0.5, electron_fraction_ocean=0.5), False),
        ('core override 0.40, r=1.0', dict(electron_fraction_core=0.40), True),
    ]
    for label, kw, expect in cases:
        got = bool(warned(**kw))
        ok = got == expect
        bad += [] if ok else [('guard: ' + label, got)]
        print('  %-38s warns=%-5s expected=%-5s %s'
              % (label, got, expect, '' if ok else '  <-- WRONG'))

    print()
    print('interactions failing:', len(bad))
    for b in bad:
        print('   ', b)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
