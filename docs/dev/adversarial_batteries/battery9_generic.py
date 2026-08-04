# -*- coding: utf-8 -*-
"""Battery 9 -- the generic, arbitrary-Hamiltonian entry points.

`_osc_prob_hybrid_dispatch_generic` serves osc_prob_sun and osc_prob_earth, which take a
user-supplied H_func(energy, l, VCC) rather than the separable vacuum/matter decomposition the
std/NSI/LIV wrappers build.  It is the same hybrid machinery reached by a different door, and
every defect the validation found lived in that machinery -- yet it had exactly one row of
coverage (a bit-identity check).  This closes that.

PASS CRITERIA, stated before running:
  (H1) every answer meets the requested tolerance OR warns.  A silent miss is a failure.
  (H2) the two safeguards added for the separable path must be inherited here, since they live
       inside hybrid_propagator: an unmarked discontinuity must not certify, and sub-threshold
       gamma must not certify on agreement alone.
  (H3) osc_prob_earth on real PREM must be correct or warn, in both directions (with the layer
       breakpoints supplied, and without).

Oracle: solve_ivp/DOP853, or expm where a profile makes expm exact.
"""

import sys
import time
import warnings

import numpy as np
from scipy.linalg import expm

import harness as H
import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.matter as matter
import magnus.oscprob as op

TOL = 1e-3
WARN_OK = ('ToleranceNotAchievedWarning', 'MagnusConvergenceWarning',
           'HybridCertificationWarning', 'UnmarkedDiscontinuityWarning')


def _hvac(d):
    p = H.params_for(d)
    return np.asarray(H.h_vac_for(d, p), dtype=complex)


def user_H(d, kind='plain'):
    """A user Hamiltonian in the shape osc_prob_sun/osc_prob_earth expect: H(energy, l, VCC).

    `kind` selects an adversarial extra term that lives in H directly rather than in the
    density -- which is the handle a user has and the wrappers do not.
    """
    hv = _hvac(d)
    proj = np.zeros((d, d), dtype=complex)
    proj[0, 0] = 1.0
    off = np.zeros((d, d), dtype=complex)
    off[0, 1] = off[1, 0] = 1.0

    def H_func(energy, l, VCC, _k=kind):
        vcc = np.asarray(VCC)
        base = (1.0/energy)*hv + vcc[..., None, None]*proj
        if _k == 'plain':
            return base
        x = np.asarray(l, dtype=float)
        scale = float(np.max(np.abs(hv)))/float(energy)
        if _k == 'step':                       # jump in H itself, nothing to do with VCC
            f = np.where(x < 0.5*x.max() if x.ndim else x < 1.0, 0.0, 1.0)
            return base + 0.30*scale*np.asarray(f)[..., None, None]*off
        if _k == 'kink':
            f = np.abs(x - 0.5*np.max(x)) if x.ndim else np.abs(x)
            return base + 0.05*scale*(np.asarray(f)/max(float(np.max(np.abs(x))), 1.0)
                                      )[..., None, None]*off
        raise ValueError(_k)

    return H_func


def step_user_H(d, mid, amp=0.30):
    """H with an unmarked jump at `mid`, independent of the density.  Constant on each side if
    the caller also passes a constant VCC, which makes expm exact there."""
    hv = _hvac(d)
    proj = np.zeros((d, d), dtype=complex)
    proj[0, 0] = 1.0
    off = np.zeros((d, d), dtype=complex)
    off[0, 1] = off[1, 0] = 1.0

    def H_func(energy, l, VCC):
        vcc = np.asarray(VCC)
        x = np.asarray(l, dtype=float)
        f = np.where(x < mid, 0.0, 1.0)
        scale = float(np.max(np.abs(hv)))/float(energy)
        return ((1.0/energy)*hv + vcc[..., None, None]*proj
                + amp*scale*np.asarray(f)[..., None, None]*off)

    return H_func


def score(label, call, H_of_l, l0, l1, d, exact=None):
    with H.Caught() as c:
        t0 = time.time()
        P = np.asarray(call())
        dt = time.time() - t0
    P = P.reshape((d, d)) if P.size == d*d else P
    P_ex = exact if exact is not None else H.P_of(H.exact_U(H_of_l, l0, l1, d))
    err = H.maxabs(P - P_ex)
    silent = err > TOL and not any(x in WARN_OK for x in c.names)
    print('  %-46s err=%9.3e  %6.2fs  %-11s %s'
          % (label, err, dt, 'SILENT MISS' if silent else ('warned' if err > TOL else 'ok'),
             ','.join(c.names) or '-'), flush=True)
    return silent


def sub1():
    """9.1 osc_prob_sun with a user Hamiltonian: smooth, and the two adversarial shapes."""
    print('## 9.1  osc_prob_sun (generic H), against solve_ivp')
    E, L = 20.0*gd.UNIT_MEV, 0.4*gd.SUN_RADIUS*gd.UNIT_KM
    vcc = H.vcc_of(H.solar_ne())
    bad = 0
    for d in (2, 3, 4):
        for kind in ('plain', 'kink'):
            Hf = user_H(d, kind)

            def one_arg(l, _H=Hf, _E=E):
                return _H(_E, l, np.asarray(vcc(l)))

            bad += score('d=%d %-6s' % (d, kind),
                         lambda _H=Hf: op.osc_prob_sun(_H, E, L, validate_input=False),
                         one_arg, 0.0, L, d)
    print()
    return bad


def sub2():
    """9.2 The two safeguards must be inherited by the generic path."""
    print('## 9.2  safeguards inherited from hybrid_propagator (criterion H2)')
    E, L = 50.0*gd.UNIT_MEV, 0.3*gd.SUN_RADIUS*gd.UNIT_KM
    mid = 0.5*L
    vcc = H.vcc_of(H.solar_ne())
    bad = 0
    for d in (2, 3):
        Hf = step_user_H(d, mid)

        def one_arg(l, _H=Hf, _E=E):
            return _H(_E, l, np.asarray(vcc(l)))

        # A jump inside H, on top of a smooth density: the density-based guard cannot see it.
        bad += score('d=%d unmarked jump INSIDE H (not in VCC)' % d,
                     lambda _H=Hf: op.osc_prob_sun(_H, E, L, validate_input=False),
                     one_arg, 0.0, L, d)
        # And the same through strategy='magnus', which should be right and warn.
        score('d=%d   same, strategy=magnus (control)' % d,
              lambda _H=Hf: op.osc_prob_sun(_H, E, L, strategy='magnus',
                                            validate_input=False),
              one_arg, 0.0, L, d)
    print()
    return bad


def sub3():
    """9.3 osc_prob_earth on real PREM -- discontinuous by nature."""
    print('## 9.3  osc_prob_earth, real PREM (criterion H3)')
    bad = 0
    for costhz in (-0.2, -0.6, -0.95):
        Lc = 2.0*gd.EARTH_RADIUS*abs(costhz)*gd.UNIT_KM
        for E_GeV in (1.0, 10.0):
            E = E_GeV*gd.UNIT_GEV
            for d in (2, 3):
                hv = _hvac(d)
                proj = np.zeros((d, d), dtype=complex)
                proj[0, 0] = 1.0

                def H_func(energy, l, VCC, _hv=hv, _p=proj):
                    return (1.0/energy)*_hv + np.asarray(VCC)[..., None, None]*_p

                # The true one-argument H, via the package's own PREM potential.
                vcc_prem = matter.vcc_func_from_rho_func(
                    lambda x, _c=costhz: earth.density_matter_func_prem(x, _c),
                    0.0, 1.0, 0.5, nubar=False, density_matter_is_in_g_per_cm3=True,
                    density_is_of_number_of_electrons=False)

                def one_arg(l, _E=E, _hv=hv, _p=proj):
                    return (1.0/_E)*_hv + np.asarray(vcc_prem(l))[..., None, None]*_p

                bad += score('costhz=%.2f E=%4.1f GeV d=%d' % (costhz, E_GeV, d),
                             lambda _H=H_func: op.osc_prob_earth(
                                 _H, E, costhz=costhz, L=Lc, validate_input=False),
                             one_arg, 0.0, Lc, d)
    print()
    return bad


def sub4():
    """9.4 A user Hamiltonian that is constant on each half: expm is EXACT, so no oracle doubt."""
    print('## 9.4  piecewise-constant user H, scored against expm (exact)')
    E = 50.0*gd.UNIT_MEV
    L = 0.2*gd.SUN_RADIUS*gd.UNIT_KM
    mid = 0.5*L
    bad = 0
    for d in (2, 3):
        hv = _hvac(d)
        off = np.zeros((d, d), dtype=complex)
        off[0, 1] = off[1, 0] = 1.0
        scale = float(np.max(np.abs(hv)))/float(E)

        def H_func(energy, l, VCC, _hv=hv, _o=off, _s=scale):
            x = np.asarray(l, dtype=float)
            f = np.where(x < mid, 0.0, 1.0)
            return (1.0/energy)*_hv + 0.4*_s*np.asarray(f)[..., None, None]*_o

        Ha = (1.0/E)*hv
        Hb = (1.0/E)*hv + 0.4*scale*off
        U = expm(-1j*Hb*(L - mid)) @ expm(-1j*Ha*mid)
        P_ex = H.P_of(U)

        # rho_func = 0 makes VCC a constant zero, so H is exactly the two-piece matrix above.
        bad += score('d=%d two-piece user H (VCC=0)' % d,
                     lambda _H=H_func: op.osc_prob_sun(_H, E, L, validate_input=False),
                     None, 0.0, L, d, exact=P_ex)
    print()
    return bad


if __name__ == '__main__':
    warnings.simplefilter('always')
    print('# Battery 9 -- generic (arbitrary-Hamiltonian) entry points.  tol=%.0e\n' % TOL)
    subs = {'1': sub1, '2': sub2, '3': sub3, '4': sub4}
    total = 0
    for w in (sys.argv[1:] or ['1', '2', '3', '4']):
        total += subs[w]()
    print('=== BATTERY 9: %d silent misses ===' % total)
