"""Feature interactions: `angles` crossed with every other axis, and the invariants
that must survive each combination.

This package's own history is the argument for doing it this way.  A 7440x accuracy hole
survived because each axis was swept alone -- clustered spectra at norm 1, and large norms
at generic separation -- and the damage needed both at once.  A newly added keyword that
behaves on its own says nothing about what it does when combined with the antineutrino
sign, the averaged limit, a different strategy, a different exponential backend, or a
cumulative scan.

Invariants checked per cell:
  * CONV  -- all four `angles` conventions give the same probability
  * UNIT  -- probabilities out of one initial flavour sum to 1
  * RANGE -- every probability is in [0, 1]
"""
import itertools
import sys
import warnings

import numpy as np

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.oscprob as op

P = gd.OSC_PARAMS_PREDEFINED['OSC_PARAMS_DEFAULT']
TH = {k: float(np.arcsin(P[k])) for k in ('s12', 's23', 's13')}
STER = {'s14': 0.15, 's15': 0.08, 's24': 0.10, 's25': 0.05, 's34': 0.05, 's35': 0.03}
STER_TH = {k: float(np.arcsin(v)) for k, v in STER.items()}
PH = {'dCP': P['dCP'], 'd14': 0.4, 'd15': 0.3, 'd24': 0.2, 'd35': 0.1}
D3 = {'D21': P['D21'], 'D31': P['D31']}
D4 = dict(D3, D41=1.0)
D5 = dict(D4, D51=2.0)

E = 1.0*gd.UNIT_GEV
L_VAC = 1300.0*gd.UNIT_KM
RHO = 2.8
COSTHZ = -0.6
L_EARTH = earth.distance_traveled_inside_earth(COSTHZ)*gd.CONV_KM_TO_INV_EV


def angles_in(conv):
    """Every angle and phase, stated in `conv`."""
    out = {}
    for k, th in list(TH.items()) + list(STER_TH.items()):
        out[k] = (float(np.sin(th)) if conv == 'sin' else
                  float(np.sin(th))**2 if conv == 'sin2' else
                  th if conv == 'rad' else float(np.degrees(th)))
    for k, v in PH.items():
        out[k] = float(np.degrees(v)) if conv == 'deg' else float(v)
    return out


def pick(d, names):
    return {n: d[n] for n in names}

A3 = ('s12', 's23', 's13', 'dCP')
A4 = A3 + ('s14', 'd14', 's24', 'd24', 's34')
A5 = A4 + ('s15', 'd15', 's25', 's35', 'd35')


def build(flavours, environment, conv, extra):
    """One call, in one convention, with `extra` carrying the other feature axes."""
    a = angles_in(conv)
    if flavours == 3:
        std, D = pick(a, A3), D3
    elif flavours == 4:
        std, D = pick(a, A4), D4
    else:
        std, D = pick(a, A5), D5
    common = dict(std, **D, angles=conv, **extra)
    if environment == 'vacuum':
        fn = getattr(op, 'osc_prob_%dnu_vacuum' % flavours)
        return fn(E, L_VAC, **common)
    if environment == 'constant':
        fn = getattr(op, 'osc_prob_%dnu_matter_constant_density' % flavours)
        return fn(E, L_VAC, RHO, density_matter_is_in_g_per_cm3=True, **common)
    fn = getattr(op, 'osc_prob_%dnu_earth' % flavours)
    return fn(E, costhz=COSTHZ, L=L_EARTH, **common)


AXES = {
    'flavours': [3, 4, 5],
    'environment': ['vacuum', 'constant', 'earth'],
    'nubar': [False, True],
    'average': [False, True],
}


def main():
    rows = list(itertools.product(*AXES.values()))
    print('%d cells x 4 conventions\n' % len(rows))
    print('%-9s %-11s %-7s %-8s %-11s %-11s %s'
          % ('flavours', 'environment', 'nubar', 'average', 'CONV', 'UNIT', 'RANGE'))
    print('-'*74)
    bad = []
    for flavours, environment, nubar, average in rows:
        extra = {'nubar': nubar, 'average': average}
        try:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                mats = {c: np.asarray(build(flavours, environment, c, extra), dtype=float)
                        for c in gd.ANGLE_CONVENTIONS}
        except Exception as exc:                       # noqa: BLE001
            print('%-9d %-11s %-7s %-8s  RAISED %s' % (flavours, environment, nubar,
                                                       average, str(exc)[:40]))
            bad.append((flavours, environment, nubar, average, 'raised'))
            continue

        ref = mats['sin']
        conv = max(float(np.max(np.abs(mats[c] - ref))) for c in gd.ANGLE_CONVENTIONS)
        unit = float(np.max(np.abs(ref.sum(axis=-1) - 1.0)))
        lo, hi = float(ref.min()), float(ref.max())
        rng = max(0.0, -lo, hi - 1.0)

        flags = []
        if conv > 1e-9:
            flags.append('CONV')
        if unit > 1e-9:
            flags.append('UNIT')
        if rng > 1e-9:
            flags.append('RANGE')
        if flags:
            bad.append((flavours, environment, nubar, average, ','.join(flags)))
        print('%-9d %-11s %-7s %-8s %-11.2e %-11.2e %-9.2e %s'
              % (flavours, environment, nubar, average, conv, unit, rng,
                 '<-- ' + ','.join(flags) if flags else ''))

    print()
    print('cells failing an invariant:', len(bad))
    for b in bad:
        print('   ', b)
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
