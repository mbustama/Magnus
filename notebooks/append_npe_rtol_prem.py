# -*- coding: utf-8 -*-
r"""append_npe_rtol_prem.py

Add a tolerance-dialled NuOscProbExact series to the PREM-chord benchmarks, beside the
slab-count one already there.

WHY NOT ``earth.probabilities_Nnu_earth``.  That is the natural call, and it does place
slabs on the PREM shells and refine to a tolerance.  But it builds its own matter
potential, and the two packages disagree there by a constant factor: measured at five
radii from 1000 to 6370 km, the PREM densities agree *exactly* (0.0e+00) while V_CC
differs by 1.90e-04 at every one of them.  That is a units convention, not physics, but
it is four orders above where these curves reach, so scoring both codes against one
reference would put a bookkeeping offset on the axis and call it method.

WHAT THIS DOES INSTEAD.  It gives NuOscProbExact the same Hamiltonian Mag(nu)s is given,
while still using NuOscProbExact's own geometry, its own composition and its own
refinement:

  * ``earth.earth_slabs(costhz, n)`` supplies the slab widths -- cut at every shell
    boundary the chord crosses, then subdivided.  This is what makes the discretisation
    converge at all; uniform slabs across a density jump do not, and the stored
    ``n_slabs`` series shows it, plateauing near 1e-4 and refusing to improve.
  * the Hamiltonian is built here from Mag(nu)s's V_CC, batched over the twelve energies,
    so both codes see one Hamiltonian and one reference serves both.
  * ``slabs.probabilities_Nnu_slabs`` composes the slabs, and ``slabs._n_for_tolerance``
    -- the same refinement loop ``probabilities_Nnu_earth`` uses internally, and
    shape-agnostic, so it accepts the batch -- turns the tolerance into a slab count.

The private ``_n_for_tolerance`` is deliberate: reimplementing its acceptance test would
make the refinement ours rather than the library's, which is the one part of this that
must not be.

WHAT IS KEPT.  The existing ``NuOscProbExact`` series, dialled by ``n_slabs``, stays.

    python notebooks/append_npe_rtol_prem.py
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import datetime
import json
import pathlib
import re
import shutil
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
sys.path.insert(0, '/home/mbustamante/Research/NuOscProb/NuOscProbExact/src')
import gen_profile_benchmarks as gpb                      # noqa: E402
import prem_chord_common as pcc                           # noqa: E402

from magnus import globaldefs as gd, matter               # noqa: E402

HERE = pathlib.Path(__file__).parent
STORE = HERE/'external_prem_chord_benchmarks.json'
N_START = 2
N_MAX = 32768          # per segment; 17 segments, so ~557k slabs and ~1.7 GB at 4nu
SERIES = 'NuOscProbExact, rtol'


def npe_rtol_points(d, E, dials, ref):
    import earth as npe_earth
    import slabs as npe
    fn = getattr(npe, 'probabilities_%dnu_slabs' % d, None)
    if fn is None:
        return []                                   # no five-flavour route exists

    hv = gpb.h_vac(d)
    proj = gpb.matter.matter_potential_projector(d)
    per_rho = matter.VCC_func(0.0, lambda l: matter.num_density_e_func(
        0.0, lambda _x: 1.0, electron_fraction=pcc.ELECTRON_FRACTION,
        density_matter_is_in_g_per_cm3=True))

    def evaluate(n):
        """Probabilities at n sub-slabs per PREM segment, batched over the energies."""
        widths_km, rho = npe_earth.earth_slabs(pcc.COSTHZ, n)
        # earth_slabs returns kilometres; probabilities_Nnu_slabs documents eV^-1.
        # Passing the raw km converges -- and is wrong: it was accepted at two slabs
        # per segment with the survival probability off by 0.92.  Magnus's constant
        # is the right one here, since it makes the widths sum to exactly the
        # baseline the reference was built on.
        widths = np.asarray(widths_km, dtype=float)*gd.CONV_KM_TO_INV_EV
        v = np.asarray(rho, dtype=float)*per_rho        # Magnus's V_CC, not NPE's
        H = np.broadcast_to((hv[None, None]/E[:, None, None, None]).astype(complex),
                            (len(E), len(v), d, d)).copy()
        H += v[None, :, None, None]*proj[None, None]
        return np.asarray(fn(H, widths), dtype=float)

    col = gd.NUMU*d + gd.NUMU
    out = []
    for rtol in dials:
        used = {}

        def call(r=rtol, used=used):
            # The library's own refinement loop, not a reimplementation of its
            # acceptance test: that is the one part of this that must be theirs.
            # atol = rtol/100 matches how gen_prem_benchmarks.py dials Magnus.
            n, p = npe._n_for_tolerance(evaluate, r, r*1.0e-2, N_START, N_MAX,
                                        'append_npe_rtol_prem')
            used['n'] = n
            return np.asarray(p)

        # One dial must never cost the others.  An earlier run lost four good
        # 2nu points because the last dial raised while parsing its own error
        # message, and the case is only written once every dial has returned.
        try:
            P = call()
        except Exception as exc:                    # noqa: BLE001 -- see above
            msg = str(exc)
            # Both message shapes end the number with a sentence period, so the
            # match must stop on a digit: '...was 2.062e-10.' parsed as a float
            # with the period attached is what killed the first run.
            best = re.search(r'was\s+([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)', msg)
            try:
                floor = float(best.group(1)) if best else None
            except ValueError:
                floor = None
            reached = isinstance(exc, ValueError) and 'could not meet' in msg
            print('    rtol %.0e  %s -- %s'
                  % (rtol, 'UNREACHABLE' if reached else type(exc).__name__,
                     ('estimate near %s' % best.group(1)) if best else msg[:70]),
                  file=sys.stderr, flush=True)
            out.append(dict(label='%.0e' % rtol, rtol=rtol, unreachable=True,
                            best_error_estimate=floor,
                            note=msg.split('.  ')[0][:200]))
            continue
        try:
            us = 1.0e6*gpb.timed(call)/len(E)
            err = float(np.max(np.abs(P[..., col] - ref)))
        except Exception as exc:                    # noqa: BLE001
            print('    rtol %.0e  timing/scoring failed: %s'
                  % (rtol, exc), file=sys.stderr, flush=True)
            continue
        out.append(dict(label='%.0e' % rtol, rtol=rtol,
                        n_slabs_per_segment=int(used['n']),
                        us_per_probability=us, max_abs_error=err))
        print('    rtol %.0e  %9.2f us/prob  err %.3e  (%d per segment)'
              % (rtol, us, err, used['n']), file=sys.stderr, flush=True)
    return out


def main():
    store = json.loads(STORE.read_text())
    shutil.copy(STORE, STORE.with_suffix('.json.pre_npe_rtol'))
    E = np.asarray(store['cases'][0]['energy_ev'], dtype=float)

    for case in store['cases']:
        d = case['flavours']
        if any(s['name'] == SERIES for s in case['series']):
            continue
        print('--- %d flavours ---' % d, file=sys.stderr, flush=True)
        pts = npe_rtol_points(d, E, gpb.MG_DIALS, np.array(case['reference']))
        if not pts:
            print('    no %dnu route' % d, file=sys.stderr, flush=True)
            continue
        case['series'].append({'name': SERIES, 'dial': 'rtol', 'points': pts})
        STORE.write_text(json.dumps(store, indent=1))

    store['npe_rtol_series'] = {
        'added': datetime.date.today().isoformat(),
        'nuoscprobexact': '1.14.0',
        'route': ('earth.earth_slabs for the shell-aligned geometry, '
                  'slabs.probabilities_Nnu_slabs to compose, '
                  'slabs._n_for_tolerance to refine'),
        'hamiltonian': ("Magnus's V_CC, shared with the Magnus series, so one reference "
                        "serves both; the packages' own V_CC differ by 1.90e-04"),
    }
    STORE.write_text(json.dumps(store, indent=1))
    print('wrote %s' % STORE.name, file=sys.stderr)


if __name__ == '__main__':
    main()
