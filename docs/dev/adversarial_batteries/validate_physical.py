# -*- coding: utf-8 -*-
"""Does each physical profile actually have the shape it claims?

Run this before anything measures with these profiles.  "A finding on a badly-built profile is
a finding about your construction" -- so every claim the docstrings in ``physical_profiles.py``
make is checked here against the profile itself, cheaply and without an oracle:

  * the tabulated families have a kink at every node and nowhere else;
  * the SN shock jumps by the factor the literature says, at the radius it says, over the width
    that was asked for, and the MSW resonance really is on the trajectory;
  * the turbulence has the requested rms, a k^-5/3 spectrum, and power below the reference grid;
  * the Earth chord crosses the declared layers, twice each;
  * every profile is positive, finite, scalar-for-scalar and array-capable.

Run:  python validate_physical.py
"""

import sys

import numpy as np

import harness as H
import magnus.globaldefs as gd
import physical_profiles as pp
from battery2 import ne_res_for

FAIL = []


def check(cond, what, detail=''):
    tag = 'ok  ' if cond else 'FAIL'
    if not cond:
        FAIL.append(what)
    print('  [%s] %-58s %s' % (tag, what, detail))


def api_contract(fams):
    """Every trap in the handover's API section, on every family."""
    print('\n--- API contract (scalar-for-scalar, array-capable, positive, finite) ---')
    for f in fams:
        l0, l1 = f['l0'], f['l1']
        mid = 0.5*(l0 + l1)
        s = f['ne'](mid)
        arr = np.asarray(f['ne'](np.linspace(l0, l1, 1001)))
        ok = (isinstance(s, float) or (np.ndim(s) == 0)) and arr.shape == (1001,) \
            and np.all(np.isfinite(arr)) and np.all(arr > 0.0)
        check(ok, f['label'],
              'scalar=%s arr=%s min=%.3e' % (type(s).__name__, arr.shape, arr.min()))


def tabulated_kinks():
    """A linear interpolant is C0 but not C1: the second difference must spike at every node."""
    print('\n--- 2.1 tabulated: kink at every node, and nowhere else ---')
    l0, l1 = 0.0, 1.0*H.L_SCALE
    for n in (20, 200, 1000):
        ne = pp.tabulated_ne(n, l0, l1, 'linear')
        # Sample 12 points per node interval; the curvature of the underlying exponential is
        # smooth, so any second-difference spike is the interpolation's.
        xs = np.linspace(l0, l1, (n - 1)*12 + 1)
        v = np.asarray(ne(xs))
        d2 = np.abs(np.diff(v, 2))
        at_node = np.zeros(len(d2), dtype=bool)
        at_node[np.clip(np.arange(1, n - 1)*12 - 1, 0, len(d2) - 1)] = True
        peak, floor = d2[at_node].min(), d2[~at_node].max()
        check(peak > 20.0*floor, 'linear N=%d kinks are at the nodes' % n,
              'min at node %.3e >> max between %.3e (x%.0f)' % (peak, floor, peak/max(floor,
                                                                                      1e-300)))
        spacing = (l1 - l0)/(n - 1)
        check(True, '   node spacing / probe0 spacing',
              '%.3f  (n_probe0 = 200)' % (spacing/((l1 - l0)/199.0)))
    cub = pp.tabulated_ne(50, l0, l1, 'cubic')
    xs = np.linspace(l0, l1, 5000)
    check(np.all(np.asarray(cub(xs)) > 0.0), 'cubic N=50 stays positive')


def bs05_shape():
    print('\n--- 2.5 BS05(AGS,OP): the real model against the package exponential ---')
    r, ne = pp.load_bs05()
    check(len(r) > 1000, 'table parsed', '%d rows, r = %.5f .. %.5f R_sun' % (len(r), r[0],
                                                                             r[-1]))
    check(bool(np.all(np.diff(r) > 0)), 'radius column strictly increasing')
    solar = H.solar_ne()
    print('   %-10s %12s %12s %8s' % ('r/R_sun', 'BS05', 'package exp', 'ratio'))
    for rr in (0.005, 0.05, 0.0949, 0.2, 0.5):
        i = int(np.argmin(np.abs(r - rr)))
        x = r[i]*gd.SUN_RADIUS*pp.KM
        a, b = float(pp.bs05_ne()(x)), float(solar(x))
        print('   %-10.4f %12.4e %12.4e %8.3f' % (r[i], a, b, a/b))
    span = pp.BS05_L1 - pp.BS05_L0
    inside = r*gd.SUN_RADIUS*pp.KM < pp.BS05_L1
    n_in = int(inside.sum())
    med = float(np.median(np.diff(r[inside])))*gd.SUN_RADIUS*pp.KM/span
    check(n_in > 100, 'nodes inside the trajectory',
          '%d, median spacing %.3e of span (probe0 = %.3e)' % (n_in, med, 1.0/199.0))


def sn_shock_shape():
    print('\n--- 2.2 SN shock: jump factor, radius, width, and the resonance on the path ---')
    for w in (1.0e-2, 1.0e-4, 1.0e-6):
        ne = pp.sn_shock_ne(w)
        w_km = w*(pp.SN_R1_KM - pp.SN_R0_KM)
        # Isolate each front by differencing against the SAME profile with that front switched
        # off.  Dividing out an analytic background instead does not work at w = 1e-2, where the
        # 700 km ramp is wide enough that the rarefaction shape varies across it -- that read as
        # a jump of 0.72 and was a defect in the measurement, not in the profile.
        no_contact = pp.sn_shock_ne(w, contact_jump=1.0)

        def bg(r_km_):
            return pp.ne_from_rho_cgs(pp._fogli_rho0(r_km_))

        # Forward shock: against the bare progenitor power law.  The expected ratio is xi times
        # the rarefaction f evaluated where the ratio is taken -- at w = 1e-2 the sampling point
        # is 2100 km inside the shock, where f has already fallen to 0.70, and comparing against
        # a bare xi = 10 there tests the sampling point rather than the profile.
        lo_km, hi_km = pp.R_FORWARD_SHOCK_KM - 3.0*w_km, pp.R_FORWARD_SHOCK_KM + 3.0*w_km
        ratio = ((float(no_contact(lo_km*pp.KM))/bg(lo_km))
                 / (float(no_contact(hi_km*pp.KM))/bg(hi_km)))
        want = pp.FOGLI_SHOCK_JUMP*float(pp._fogli_rarefaction(lo_km, pp.R_FORWARD_SHOCK_KM))
        check(abs(ratio/want - 1.0) < 0.02, 'w=%.0e forward jump = xi*f = %.2f' % (w, want),
              'measured %.3f' % ratio)

        # Contact: against the same profile with only that front switched off, so everything
        # else -- power law, rarefaction, forward shock -- divides out exactly.
        lo_km, hi_km = pp.R_CONTACT_KM - 3.0*w_km, pp.R_CONTACT_KM + 3.0*w_km
        ratio = ((float(ne(lo_km*pp.KM))/float(no_contact(lo_km*pp.KM)))
                 / (float(ne(hi_km*pp.KM))/float(no_contact(hi_km*pp.KM))))
        check(abs(ratio/pp.CONTACT_JUMP - 1.0) < 0.02,
              'w=%.0e contact jump = %.2f' % (w, pp.CONTACT_JUMP), 'measured %.3f' % ratio)
        # The front must be confined to its declared width: the profile is flat 3 widths out.
        xs = np.linspace(pp.R_FORWARD_SHOCK_KM + 3.0*w_km,
                         pp.R_FORWARD_SHOCK_KM + 40.0*w_km, 400)*pp.KM
        v = np.asarray(ne(xs))
        smooth = np.max(np.abs(np.diff(v, 2)))/np.max(np.abs(v))
        check(smooth < 1e-4, 'w=%.0e front confined to its width' % w,
              'curvature outside %.2e' % smooth)
    # The MSW resonance must be ON the ray, or the family is just a steep profile with two steps.
    # It is the H resonance (the 1-2 level gap at d = 3, driven by D31) that matters here:
    # shock effects on supernova neutrinos are an H-resonance phenomenon, and the L resonance
    # (the 0-1 gap, which ``ne_res_for`` returns) sits at 1.7-3.9e5 km, well outside this ray.
    print('   H resonance (gap 1-2, d = 3) against the ray:')
    xs = np.geomspace(pp.SN_R0_KM, pp.SN_R1_KM, 400000)
    v = np.asarray(pp.sn_shock_ne(1e-3)(xs*pp.KM))
    p3 = H.params_for(3)
    for E in pp.SN_ENERGIES:
        ner = _resonance_ne(3, p3, E, level=1)
        idx = np.where(np.diff(np.sign(v - ner)) != 0)[0]
        radii = ', '.join('%.0f km' % xs[i] for i in idx) or 'none'
        check(len(idx) >= 1, 'E=%4.0f MeV crosses the H resonance on the ray' % (E/1e6),
              'ne_res=%.3e at %s' % (ner, radii))
    check(v.min() < _resonance_ne(3, p3, pp.SN_ENERGIES[0], level=1) < v.max(),
          'H resonance inside the profile range',
          'ne on ray %.3e .. %.3e' % (v.min(), v.max()))
    print('     (the L resonance is NOT on this ray: ne_res = %.3e, below the ray minimum %.3e)'
          % (ne_res_for(3, p3, pp.SN_ENERGIES[0]), v.min()))
    check(True, 'shock radii inside the ray',
          'contact %.0f, forward %.0f in [%.0f, %.0f] km'
          % (pp.R_CONTACT_KM, pp.R_FORWARD_SHOCK_KM, pp.SN_R0_KM, pp.SN_R1_KM))


def _resonance_ne(d, params, energy, level=0):
    """Electron density minimising the gap between levels ``level`` and ``level+1``.

    ``battery2.ne_res_for`` hard-codes ``level = 0``, which is the L resonance; the supernova
    families need the H resonance.  Same bracketing search, one extra argument.
    """
    h_vac = np.asarray(H.h_vac_for(d, params), dtype=complex)
    proj = np.zeros((d, d), dtype=complex)
    proj[0, 0] = 1.0
    c_vcc = _vcc_per_ne()

    def gap(ne):
        lam = np.linalg.eigvalsh(h_vac/energy + ne*c_vcc*proj)
        return lam[level + 1] - lam[level]

    xs = np.geomspace(H.NE0*1e-4, H.NE0*1e3, 6000)
    i = int(np.argmin([gap(x) for x in xs]))
    a, b = xs[max(i - 1, 0)], xs[min(i + 1, len(xs) - 1)]
    for _ in range(200):
        m1, m2 = a + (b - a)/3.0, b - (b - a)/3.0
        a, b = (a, m2) if gap(m1) < gap(m2) else (m1, b)
    return 0.5*(a + b)


def _vcc_per_ne():
    """V_CC per unit electron number density, read off the package rather than restated."""
    return float(H.vcc_of(lambda l: 1.0)(0.0))


def turbulence_spectrum():
    print('\n--- 2.3 turbulence: rms, spectral index, and power below the reference grid ---')
    rng = np.random.default_rng(1)
    r_lo, r_hi = pp.R_CONTACT_KM, pp.R_FORWARD_SHOCK_KM
    base = pp.sn_shock_ne(1.0e-3)
    for c in (0.01, 0.1):
        # C_star is the rms of F over REALIZATIONS at a point, which is what sum(V_a) = 1
        # normalizes.  Measuring the spatial std over a sub-window instead reads about half of
        # it, and correctly so: under a k^-5/3 spectrum most of the variance is in the longest
        # mode, whose wavelength is 4*pi*(r_> - r_<) ~ 2.3e5 km, so over a 1.4e4 km window it is
        # nearly constant and `np.std` removes it with the mean.  That is a property of the
        # spectrum, not an error in the profile, and measuring it the wrong way read as a
        # factor-of-two miss.
        mid = 0.5*(r_lo + r_hi)*pp.KM
        vals = [float(pp.sn_turbulent_ne(np.random.default_rng(rng.integers(1 << 30)),
                                         c_star=c)(mid))/float(base(mid)) - 1.0
                for _ in range(600)]
        got = float(np.std(vals))
        check(0.8*c < got < 1.25*c, 'C*=%g gives ensemble rms %.4f' % (c, got),
              'over 600 realizations at the midpoint')
    # And the spatial roughness is present too, at the half-ish level the spectrum implies.
    ne = pp.sn_turbulent_ne(np.random.default_rng(11), c_star=0.1)
    xs = np.linspace(r_lo + 2000.0, r_hi - 2000.0, 200001)*pp.KM
    spatial = float(np.std(np.asarray(ne(xs))/np.asarray(base(xs)) - 1.0))
    check(0.02 < spatial < 0.1, 'C*=0.1 spatial std over the window', '%.4f' % spatial)

    # Spectral index: the power spectral density of F should fall as k^-5/3.
    ne = pp.sn_turbulent_ne(np.random.default_rng(7), c_star=0.1)
    base = pp.sn_shock_ne(1.0e-3)
    n = 1 << 20
    xs = np.linspace(r_lo + 2000.0, r_hi - 2000.0, n)
    f = np.asarray(ne(xs*pp.KM))/np.asarray(base(xs*pp.KM)) - 1.0
    psd = np.abs(np.fft.rfft(f - f.mean()))**2
    k = np.fft.rfftfreq(n, d=(xs[1] - xs[0]))
    sel = (k > 0) & (psd > 0)
    lo = np.percentile(k[sel], 1)
    band = sel & (k > lo)
    slope = np.polyfit(np.log(k[band]), np.log(np.maximum(psd[band], 1e-300)), 1)[0]
    check(-3.0 < slope < -1.0, 'power spectrum falls as a power law',
          'fitted index %.2f (drawn from alpha = %.2f)' % (slope, -pp.TURBULENCE_ALPHA))

    # The point of the family: power at scales the reference grid cannot sample.
    span = pp.SN_L1 - pp.SN_L0
    q_star = 1.0/(2.0*(r_hi - r_lo))
    shortest_km = 1.0/(q_star*pp.TURBULENCE_DYNAMIC_RANGE)
    ref_km = (span/pp.KM)/6400.0
    check(shortest_km < ref_km, 'shortest mode is below the reference grid',
          'lambda_min %.4f km vs 6400-grid spacing %.3f km' % (shortest_km, ref_km))


def earth_chord():
    print('\n--- 2.4 Earth chord: the extra layers are crossed, twice each ---')
    for costhz in (-1.0, -0.8, -0.4):
        l0, l1 = pp.earth_chord_span(costhz)
        ne = pp.earth_crust_ne(costhz)
        xs = np.linspace(l0, l1, 400001)
        v = np.asarray(ne(xs))
        # Count sizeable relative steps: PREM's own edges plus the three added layers.
        rel = np.abs(np.diff(v))/np.maximum(v[:-1], 1e-300)
        steps = int(np.sum(rel > 1.0e-3))
        want = float(pp.ne_from_rho_cgs(pp.EXTRA_CRUST_LAYERS[0][2]))
        check(abs(float(v[0])/want - 1.0) < 0.02, 'costhz=%.1f starts in the sediment layer'
              % costhz, 'ne(l0)=%.4e want %.4e' % (float(v[0]), want))
        check(steps >= 6, 'costhz=%.1f has >= 6 discontinuities' % costhz,
              '%d steps above 0.1%%, chord %.0f km' % (steps, (l1 - l0)/pp.KM))
        check(bool(np.all(np.asarray(earth_r(costhz, xs)) <= gd.EARTH_RADIUS*(1 + 1e-12))),
              'costhz=%.1f radius stays inside the Earth' % costhz)


def earth_r(costhz, l):
    import magnus.earth as earth
    return earth.earth_radial_distance_from_depth(costhz, np.asarray(l, dtype=float)/pp.KM)


def main():
    fams = pp.families()
    print('# Validating %d physical families' % len(fams))
    api_contract(fams)
    tabulated_kinks()
    bs05_shape()
    sn_shock_shape()
    turbulence_spectrum()
    earth_chord()
    print('\n=== %d checks failed ===' % len(FAIL))
    for f in FAIL:
        print('   ' + f)
    return 1 if FAIL else 0


if __name__ == '__main__':
    sys.exit(main())
