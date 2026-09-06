# -*- coding: utf-8 -*-
r"""gen_solar_average_cost.py

What the *averaged* solar probability costs, per physics configuration, on the
real BS2005-AGS,OP model.  Writes ``external_solar_average_cost.json``.

    python notebooks/gen_solar_average_cost.py > notebooks/external_solar_average_cost.json

WHY COST AND NOT SPEED-AGAINST-ACCURACY.  On this profile the averaged
probability has no accuracy dial.  Every eigenvalue pair is fully decohered and
the evolution is adiabatic throughout -- the adiabaticity parameter stays a
factor of eleven or more under the detector's threshold at every flavour count
from 0.1 to 20 MeV -- so the crossing matrix is the identity and the answer
reduces to two eigendecompositions, exact to ~1e-15.  Sweeping ``n_points``,
``threshold``, ``n_probe``, ``magnus_exp_order`` or ``integration_method``
returns bit-identical output.  There is therefore nothing to put on an accuracy
axis, and what remains worth showing is the price of the answer across the
scenarios the package covers.  See ``docs/dev/SCOPE_FIG13.md``.

THE PROFILE IS THE TABULATED MODEL, NOT THE EXPONENTIAL FIT.  ``osc_prob_*_sun``
uses the fit; this script uses ``gen_profile_benchmarks.solar_profile()``, which
log-interpolates the BS2005-AGS,OP table.  The fit is high by a factor 2.4 inside
0.05 R_sun, so they are different problems rather than two versions of one -- and
they cost differently: 5-8 ms on the fit against 8-40 ms here.

Timings are best-of-N wall clock with the first call discarded (it carries the
numba compile), and an interleaved control workload is recorded so a later run
can say whether this machine still times the same way.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import json
import pathlib
import platform
import sys
import time
import warnings

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import gen_profile_benchmarks as gpb                      # noqa: E402
import magnus.globaldefs as gd                            # noqa: E402
import magnus.hamiltonians as hams                        # noqa: E402
import magnus.oscprob as oscprob                          # noqa: E402

ENERGY = 5.0e6                                            # 5 MeV, the solar workhorse
REF_MAX_SECONDS = 1800.0     # anything projecting past half an hour is not run
REPEAT = 7

# Large enough to matter, small enough to stay a perturbation -- the same scale
# gen_shock_nsi.py uses, so the two figures' NSI mean the same thing.
_EPS_ACTIVE = dict(eps_ee=0.15, eps_em=0.05, eps_et=0.0,
                   eps_mm=0.0, eps_mt=0.0, eps_tt=0.0)
# The sterile sector carries no NSI of its own: these rows ask what the *same*
# active-sector NSI costs once there are more states to carry it, not what a
# richer NSI model costs.  Every added entry is therefore zero, and the names
# differ by flavour count because the wrappers unpack them positionally.
_EPS_STERILE = {4: ('eps_es', 'eps_ms', 'eps_ts', 'eps_ss'),
                5: ('eps_es1', 'eps_es2', 'eps_ms1', 'eps_ms2', 'eps_ts1',
                    'eps_ts2', 'eps_s1s1', 'eps_s1s2', 'eps_s2s2')}


def nsi_params_for(d):
    eps = dict(_EPS_ACTIVE)
    eps.update({k: 0.0 for k in _EPS_STERILE.get(d, ())})
    return eps
LIV = dict(sxi12=0.1, sxi23=0.1, sxi13=0.0, dxiCP=0.0,
           b1=gd.B1, b2=gd.B2, b3=gd.B3, Lambda=gd.LAMBDA, n_liv=1)
STERILE_KEYS = ('s14', 's24', 's34', 'd14', 'd24', 'D41')


def timed(call):
    """Best of REPEAT, first pass discarded: it carries the one-off numba compile."""
    call()
    best = np.inf
    for _ in range(REPEAT):
        t0 = time.perf_counter()
        call()
        best = min(best, time.perf_counter() - t0)
    return best


NSI_H = {3: hams.hamiltonian_3nu_nsi, 4: hams.hamiltonian_4nu_nsi,
         5: hams.hamiltonian_5nu_nsi}
LIV_H = {3: hams.hamiltonian_3nu_liv}


def h_of_case(kind, d, vcc, energy):
    """H(l) for one case, matching what that case's ``osc_prob_*`` call integrates.

    Each of these was checked against its own ``osc_prob_*`` call on a short segment
    before any long solve was launched, and the NSI arm failed that check the first
    time: ``hamiltonian_*nu_nsi`` returns the epsilon deviation ALONE -- at eps = 0 it
    is the zero matrix, not ``V*diag(1,0,0)`` -- so the standard matter term has to be
    added beside it.  Without that the reference neutrino crosses the Sun seeing only
    new physics and no ordinary matter, which agreed with nothing and would have cost
    an hour and a half of solves to discover afterwards.
    """
    hv = np.asarray(gpb.h_vac(d))/energy
    proj = gpb.matter.matter_potential_projector(d)
    if kind == 'liv':
        hl = np.asarray(LIV_H[d](energy, **LIV))
        return lambda x: hv + float(vcc(x))*proj + hl
    if kind == 'nsi':
        eps = nsi_params_for(d)
        return lambda x: (hv + float(vcc(x))*proj
                          + np.asarray(NSI_H[d](float(vcc(x)), **eps)))
    return lambda x: hv + float(vcc(x))*proj


def project_reference_seconds(hf, d, L, frac=4.0e-5):
    """What the full reference solve would cost, from a short segment of the same ray.

    The sterile rows cannot be run: the solver has to resolve the fastest phase in H,
    and Delta m^2_41 = 1 eV^2 is some 400 times the atmospheric splitting.  On an
    identical segment 3+1 takes 305162 right-hand sides against 926 at three flavours,
    which is a 330-fold step count on the same path, and the full solve runs for days.

    So the price is projected instead, and the projection is measured rather than
    argued.  What licenses it is that the step count grows LINEARLY with path length --
    doubling the segment doubles it, 1.97 at four flavours and 1.98 at five -- and that
    the same extrapolation reproduces the rows that were run: it predicts 354 s at three
    flavours against 354 s measured, and 188 s against 181 s at two.  The per-step cost
    is flat at about 16 microseconds from two flavours to five, so none of this is a
    costlier Hamiltonian at larger d; it is step count alone.

    Checked before trusting: an earlier projection scaled the eigenvalue gap of H
    instead, and while it agreed to within about 20 percent it was never scored against
    a measured row.  This one is.
    """
    def rhs(x, y):
        U = y.reshape(d, d, 2)
        dU = -1j*np.asarray(hf(x)) @ (U[..., 0] + 1j*U[..., 1])
        return np.stack([dU.real, dU.imag], axis=-1).ravel()

    y0 = np.stack([np.eye(d), np.zeros((d, d))], axis=-1).ravel()
    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (0.0, L*frac), y0, method='DOP853', rtol=1.0e-8, atol=1.0e-10)
    el = time.perf_counter() - t0
    return el/frac, int(sol.nfev)/frac


def coherent_reference_seconds(hf, d, L, rtol=1.0e-8, frac=0.01, n_avg=20001):
    """Wall clock for the AVERAGED probability by direct integration.

    This times *the same quantity the closed form returns*, not a step toward it.
    One DOP853 solve with dense output carries the coherent evolution along the ray;
    the phase average is then taken from the interpolant over a window at the far
    end, wide enough that every oscillation wraps many times inside it.  Both the
    solve and the averaging are inside the clock, because both are needed to get the
    number.

    An earlier version timed a single solve and called the result a lower bound.  It
    was one, but a lower bound on the alternative is not the alternative, and putting
    it on the same axis as the closed form invited the reader to compare two
    different quantities.  Dense output makes the honest comparison affordable: the
    window costs interpolant evaluations, not further integrations.

    The tolerance is deliberately loose.  The averaged answer is exact to about
    1e-15, so pricing the coherent route against *that* would charge it for a
    precision nobody asking for an average would request.  1e-4 is the most generous
    setting the comparison can be given.  The integrand is the same Hamiltonian the
    averaged call uses, so the two prices are of one physical quantity on one
    profile.

    Returns the seconds, the rhs count, whether the solve converged, and the averaged
    P_ee -- that last so the caller can check this arm against the closed form.  The
    cost claim only means anything if the two agree on the answer.
    """
    def rhs(x, y):
        U = y.reshape(d, d, 2)
        dU = -1j*np.asarray(hf(x)) @ (U[..., 0] + 1j*U[..., 1])
        return np.stack([dU.real, dU.imag], axis=-1).ravel()

    y0 = np.stack([np.eye(d), np.zeros((d, d))], axis=-1).ravel()
    # The window sits at the far end of the ray and spans `frac` of it.  Its sample
    # points are known before the solve, so they go in as `t_eval` rather than through
    # `dense_output`: the solver fills them as it passes and throws each step's
    # interpolant away.  `dense_output=True` keeps every one of them instead, and
    # across the millions of steps these solves take that reached 9.5 GB and was
    # OOM-killed on the four-flavour row.  Same samples, bounded memory.
    xs = np.linspace(L*(1.0 - frac), L, n_avg)
    t0 = time.perf_counter()
    sol = solve_ivp(rhs, (0.0, L), y0, method='DOP853', rtol=rtol,
                    atol=rtol*1.0e-2, t_eval=xs)
    ys = sol.y.reshape(d, d, 2, -1)
    U = ys[..., 0, :] + 1j*ys[..., 1, :]
    # A plain mean over a window holding N cycles of the SLOWEST phase leaves an
    # O(1/N) endpoint residual.  At three flavours the solar splitting gives only
    # some tens of cycles here and that residual was 2.3e-3 -- forty times the two-
    # flavour figure, and far too coarse to call the two routes agreed.  A Hann taper
    # suppresses the endpoint term at the same width and the same cost; the window
    # cannot simply be widened, since that walks the endpoint back into denser plasma
    # and quietly changes the problem being priced.
    w = np.hanning(U.shape[-1])
    p_ee = float(np.sum(w*np.abs(U[0, 0, :])**2)/np.sum(w))
    seconds = time.perf_counter() - t0
    return seconds, int(sol.nfev), bool(sol.success), p_ee


def control():
    a = np.random.default_rng(0).normal(size=(180, 180))
    return a @ a


def main():
    warnings.simplefilter('ignore')
    prof = gpb.solar_profile()
    per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)
    L = prof['baseline']

    def ne_of(x):
        return prof['vcc'](x)/per_ne

    common = dict(L0=0.0, density_is_of_number_of_electrons=True,
                  average=True, validate_input=False)

    def std(d):
        return lambda: oscprob.osc_prob_matter_std_potential(
            d, ne_of, ENERGY, L, gpb.osc_params(d), **common)

    def nsi(d):
        return lambda: oscprob.osc_prob_matter_nsi(
            d, ne_of, ENERGY, L, gpb.osc_params(d), nsi_params_for(d), **common)

    def liv(d):
        # rho_func is the sixth positional argument; omitting it silently gives
        # vacuum LIV, which runs a hundred times faster and is a different case.
        return lambda: oscprob.osc_prob_liv(
            d, ENERGY, L, gpb.osc_params(d), LIV, ne_of, **common)

    # The fourth field is the Hamiltonian family the reference solve has to build.
    # An earlier version measured one reference per flavour count and shared it across
    # that count's new-physics rows, on the assumption that NSI and LIV move the matrix
    # entries without changing how many cycles have to be resolved.  That assumption is
    # not tested anywhere, and the figure's whole claim is cost *per configuration*, so
    # every row now pays for its own solve.
    cases = [('2nu', 2, 'std', std(2)), ('3nu', 3, 'std', std(3)),
             ('3+1', 4, 'std', std(4)), ('3+2', 5, 'std', std(5)),
             ('3nu + NSI', 3, 'nsi', nsi(3)), ('3nu + LIV', 3, 'liv', liv(3)),
             ('3+1 + NSI', 4, 'nsi', nsi(4)), ('3+2 + NSI', 5, 'nsi', nsi(5))]

    REF_SECONDS = {}
    best = {'a': np.inf, 'b': np.inf}
    for _ in range(9):
        for k in ('a', 'b'):
            best[k] = min(best[k], timed(control))
    out = {'note': 'Cost of the averaged probability, BS2005-AGS,OP, per configuration',
           'profile': 'BS2005-AGS,OP, log-interpolated (NOT the exponential fit)',
           'energy_ev': ENERGY, 'baseline_inv_ev': L,
           'why_no_accuracy_axis': ('adiabatic throughout and fully decohered, so the '
                                    'crossing matrix is the identity and no knob moves '
                                    'the value; see docs/dev/SCOPE_FIG13.md'),
           'machine': platform.platform(), 'python': platform.python_version(),
           'numpy': np.__version__, 'magnus': getattr(__import__('magnus'), '__version__', '?'),
           'control_ratio': best['a']/best['b'], 'repeat': REPEAT,
           'reference_note': ('reference_seconds is the SAME averaged probability by '
                              'direct integration: one DOP853 solve at rtol 1e-8 with '
                              'dense output, then the phase average taken from the '
                              'interpolant over the last 1 percent of the ray.  Both '
                              'the solve and the averaging are inside the clock, and every '
                              'row pays for its own solve on its own Hamiltonian.  '
                              'reference_p_ee is that arm answer; compare it with p_ee '
                              'to see the two routes agree before believing the cost'),
           'cases': []}

    # Checkpoint after every row.  A full pass is about ninety minutes and the whole
    # of it used to live in memory until the final json.dump, so a session teardown
    # threw away eleven minutes of solves that had already succeeded.  Rows already
    # present in the checkpoint are not recomputed, which also makes it cheap to add
    # one configuration later without repaying for the other seven.
    ckpt = (pathlib.Path(__file__).resolve().parent
            / 'external_solar_average_cost.partial.json')
    done = {}
    if ckpt.exists():
        prev = json.loads(ckpt.read_text())
        done = {r['label']: r for r in prev.get('cases', [])}
        if done:
            print('  resuming, %d row(s) already measured: %s'
                  % (len(done), ', '.join(done)), file=sys.stderr, flush=True)

    for label, d, kind, call in cases:
        if label in done:
            out['cases'].append(done[label])
            print('  %-11s (from checkpoint)' % label, file=sys.stderr, flush=True)
            continue
        try:
            P = np.asarray(call())
        except Exception as exc:                          # noqa: BLE001
            print('  %-11s SKIPPED: %s' % (label, exc), file=sys.stderr, flush=True)
            continue
        t = timed(call)
        row = dict(label=label, flavors=d, seconds=t, ms=1.0e3*t,
                   p_ee=float(np.asarray(P)[0, 0]),
                   unitarity=float(np.max(np.abs(np.asarray(P).sum(axis=1) - 1.0))))
        # The coherent route, on this row's own Hamiltonian -- where it can be run at
        # all.  Rows whose projected solve exceeds REF_MAX_SECONDS are not attempted;
        # the sterile ones project to 39 and 78 hours apiece, and a figure cannot wait
        # on 230 hours of integration.  Those rows carry a projection instead, scaled
        # from the rows that WERE measured in this same run, and are flagged so the
        # plot can draw them as bounds rather than as data.
        hf = h_of_case(kind, d, prof['vcc'], ENERGY)
        proj, proj_nfev = project_reference_seconds(hf, d, L)
        row['reference_projected_seconds'] = proj
        row['reference_projected_nfev'] = proj_nfev
        if label not in REF_SECONDS and proj <= REF_MAX_SECONDS:
            rs, nfev, ok, ref_pee = coherent_reference_seconds(hf, d, L)
            REF_SECONDS[label] = dict(seconds=rs, nfev=nfev, converged=ok, p_ee=ref_pee)
            print('      coherent DOP853 rtol 1e-8, %-11s %8.1f s (%d rhs)  P_ee = %.6f'
                  % (label + ':', rs, nfev, ref_pee), file=sys.stderr, flush=True)
        if label in REF_SECONDS:
            row['reference_measured'] = True
            row['reference_seconds'] = REF_SECONDS[label]['seconds']
            row['reference_nfev'] = REF_SECONDS[label]['nfev']
            row['reference_p_ee'] = REF_SECONDS[label]['p_ee']
            # The cost comparison is only honest if both arms landed on the same number.
            row['reference_agreement'] = abs(row['reference_p_ee'] - row['p_ee'])
        else:
            row['reference_measured'] = False
            row['reference_seconds'] = proj
            print('      %-11s reference NOT RUN: projects to %.1f h (%.3g rhs)'
                  % (label + ':', proj/3600.0, proj_nfev), file=sys.stderr, flush=True)
        out['cases'].append(row)
        ckpt.write_text(json.dumps(out, indent=1))
        if row['reference_measured']:
            tail = 'ref %8.1f s  agree %.1e' % (row['reference_seconds'],
                                                row['reference_agreement'])
        else:
            tail = 'ref %8.1f h  PROJECTED' % (row['reference_seconds']/3600.0)
        print('  %-11s %8.2f ms   P_ee = %.6f   unitarity %.1e   %s'
              % (label, row['ms'], row['p_ee'], row['unitarity'], tail),
              file=sys.stderr, flush=True)

    json.dump(out, sys.stdout, indent=1)


if __name__ == '__main__':
    main()
