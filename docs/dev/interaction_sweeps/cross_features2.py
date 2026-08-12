"""Second interaction sweep: strategy, exponential backend, and the antineutrino identity.

Each of these is a route that must not change the answer.  `strategy` picks a propagator,
`expm_backend` picks how a 2x2/3x3 exponential is formed, and `nubar` is CP conjugation --
in vacuum that is exactly dCP -> -dCP, so the wrapper and the hand-conjugated parameter set
have to agree.  Crossed with the flavour count and the environment, because the sterile
sector and the matter term are where this package's failures have lived.
"""
import itertools
import sys
import warnings

import numpy as np

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.oscprob as op

P = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
D3 = {'D21': P['D21'], 'D31': P['D31']}
D4 = dict(D3, D41=1.0)
STD = {'s12': P['s12'], 's23': P['s23'], 's13': P['s13'], 'dCP': P['dCP']}
STER = {'s14': 0.15, 'd14': 0.0, 's24': 0.10, 'd24': 0.0, 's34': 0.05}
E = np.logspace(0.0, 1.0, 12)*gd.UNIT_GEV
L_VAC = 1300.0*gd.UNIT_KM
RHO = 2.8
COSTHZ = -0.7
L_EARTH = earth.distance_traveled_inside_earth(COSTHZ)*gd.CONV_KM_TO_INV_EV
# matched so the sterile projector and the density describe one medium, keeping this
# sweep about the routes rather than about the composition warning
R_MATCHED = 1.0


def call(flavours, environment, **extra):
    kw = dict(STD, **(D3 if flavours == 3 else D4))
    if flavours == 4:
        kw.update(STER)
        kw['electron_fraction'] = 0.5
        kw['ratio_number_neutrons_to_protons'] = R_MATCHED
    if environment == 'vacuum':
        kw.pop('electron_fraction', None)
        kw.pop('ratio_number_neutrons_to_protons', None)
        fn = getattr(op, 'osc_prob_%dnu_vacuum' % flavours)
        return np.asarray(fn(E, L_VAC, **kw, **extra), dtype=float)
    if environment == 'constant':
        fn = getattr(op, 'osc_prob_%dnu_matter_constant_density' % flavours)
        return np.asarray(fn(E, L_VAC, RHO, density_matter_is_in_g_per_cm3=True,
                             **kw, **extra), dtype=float)
    fn = getattr(op, 'osc_prob_%dnu_earth' % flavours)
    return np.asarray(fn(E, costhz=COSTHZ, L=L_EARTH, **kw, **extra), dtype=float)


def worst(a, b):
    return float(np.max(np.abs(a - b)))


def main():
    bad = []
    print('=== strategy: auto / hybrid / magnus must agree ===')
    # vacuum has no matter potential, so no adiabatic route and no `strategy` -- the
    # wrapper says so explicitly, which is why it is excluded rather than tolerated
    for flavours, env in itertools.product((3, 4), ('constant', 'earth')):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ref = call(flavours, env, strategy='magnus')
            got = {s: call(flavours, env, strategy=s) for s in ('auto', 'hybrid')}
        w = max(worst(v, ref) for v in got.values())
        flag = '' if w < 1e-6 else '  <-- DISAGREE'
        if flag:
            bad.append(('strategy', flavours, env, w))
        print('  %dnu %-9s worst |dP| vs magnus: %.3e%s' % (flavours, env, w, flag))

    print()
    print('=== expm_backend: the compiled kernel against eigh ===')
    for flavours, env in itertools.product((3, 4), ('vacuum', 'constant', 'earth')):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            a = call(flavours, env, expm_backend='eigh')
            b = call(flavours, env, expm_backend='auto')
        w = worst(a, b)
        flag = '' if w < 1e-9 else '  <-- DISAGREE'
        if flag:
            bad.append(('backend', flavours, env, w))
        print('  %dnu %-9s worst |dP| eigh vs auto: %.3e%s' % (flavours, env, w, flag))

    print()
    print('=== nubar in vacuum is exactly dCP -> -dCP ===')
    for flavours in (3, 4):
        kw = dict(STD, **(D3 if flavours == 3 else D4))
        if flavours == 4:
            kw.update(STER)
        fn = getattr(op, 'osc_prob_%dnu_vacuum' % flavours)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            a = np.asarray(fn(E, L_VAC, **kw, nubar=True), dtype=float)
            flipped = dict(kw, dCP=-kw['dCP'])
            if flavours == 4:
                flipped.update(d14=-kw['d14'], d24=-kw['d24'])
            b = np.asarray(fn(E, L_VAC, **flipped), dtype=float)
        w = worst(a, b)
        flag = '' if w < 1e-12 else '  <-- DISAGREE'
        if flag:
            bad.append(('nubar', flavours, 'vacuum', w))
        print('  %dnu vacuum   worst |P(nubar) - P(-dCP)|: %.3e%s' % (flavours, w, flag))

    print()
    print('=== integration_method: gl / trapezoid / simpson at tight tolerance ===')
    for flavours in (3,):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ref = call(flavours, 'earth', integration_method='gl', rtol=1e-8, atol=1e-10)
            for m in ('trapezoid', 'simpson'):
                got = call(flavours, 'earth', integration_method=m, rtol=1e-8, atol=1e-10)
                w = worst(got, ref)
                flag = '' if w < 1e-5 else '  <-- DISAGREE'
                if flag:
                    bad.append(('integration', flavours, m, w))
                print('  %dnu earth    %-10s vs gl: %.3e%s' % (flavours, m, w, flag))

    print()
    print('interactions failing:', len(bad))
    for b in bad:
        print('   ', b)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
