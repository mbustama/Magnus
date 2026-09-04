# -*- coding: utf-8 -*-
r"""gen_prem_benchmarks.py

Speed against accuracy on a PREM chord at cos(theta_z) = -0.9: the Earth analogue of the
smooth-profile comparison behind Figure 11.

Same twelve energies, same dials, same timing protocol and same call as that figure, so
the two are read against each other.  What differs is the profile: exponential and smooth
there, PREM and piecewise-smooth here, with sixteen layer discontinuities along the chord.
That is the comparison -- what a jump costs each method.

THE REFEREE is ``prem_chord_reference.json``: a segment-aligned triple-Richardson mpmath
ladder, cross-checked against an independent DOP853 integration.  Not the DOP853 run
alone, which on the smooth profile turned out to *be* the floor every series settled onto.

THE SIXTEEN LAYER CROSSINGS ARE DECLARED, via ``t_breakpoints``.  That is a different
keyword from ``t_slab_edges``: the latter fixes the discretization outright and would
leave rtol nothing to do, while the former only forces those positions to be slab
boundaries and lets the adaptive refinement subdivide between them.  The Earth entry
points do exactly this (``oscprob.py``, ``t_breakpoints`` from
``prem_layer_edges_along_chord``), for the reason Section 4.5 gives: the quadrature
reaches its nominal order only if the Hamiltonian is smooth inside each slab.

Measured, withholding them is not a small effect.  At 3nu and 2 GeV, rtol 1e-8 without
breakpoints returns an error of 9.8e-06 in about 18 s, having exhausted its slab cap;
with them it returns 5.2e-11 in 0.01 s.  A thousand times faster and five orders more
accurate, so a figure built without them would measure a documented failure mode rather
than the method.

    python notebooks/gen_prem_benchmarks.py

Writes ``external_prem_chord_benchmarks.json``.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"

import datetime
import json
import pathlib
import platform
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).parent))
import gen_profile_benchmarks as gpb                      # noqa: E402
import prem_chord_common as pcc                           # noqa: E402

from magnus import globaldefs as gd                       # noqa: E402

HERE = pathlib.Path(__file__).parent
REF = HERE/'prem_chord_reference.json'
OUT = HERE/'external_prem_chord_benchmarks.json'
ORDERS = (4, 6, 8)


def magnus_points(d, prof, dials, ref, order, breakpoints):
    per_ne = gpb.matter.VCC_func(l=0.0, num_density_e_func=lambda l: 1.0)

    def ne_of(x):
        return prof['vcc'](x)/per_ne

    E = prof['energies']
    osc = gpb.osc_params(d)
    out = []
    for rtol in dials:
        def call(r=rtol):
            return np.asarray(gpb.oscprob.osc_prob_matter_std_potential(
                d, ne_of, E, prof['baseline'], osc, L0=0.0,
                nu_i=gd.NUMU, nu_f=gd.NUMU,
                density_is_of_number_of_electrons=True, rtol=r, atol=r*1.0e-2,
                magnus_exp_order=order, strategy='magnus',
                t_breakpoints=breakpoints))

        P = call()
        us = 1.0e6*gpb.timed(call)/len(E)
        err = float(np.max(np.abs(P - ref)))
        out.append(dict(label='%.0e' % rtol, rtol=rtol,
                        us_per_probability=us, max_abs_error=err))
        print('      rtol %.0e  %10.2f us/prob  err %.3e' % (rtol, us, err),
              file=sys.stderr, flush=True)
    return out


def main():
    if not REF.exists():
        raise SystemExit('run gen_prem_reference.py first: %s is missing' % REF.name)
    refs = json.loads(REF.read_text())
    ch = pcc.chord()
    prof = {'vcc': ch['vcc'], 'baseline': ch['baseline'],
            'energies': np.asarray(refs['energy_ev'], dtype=float)}

    best = {'a': np.inf, 'b': np.inf}
    for _ in range(9):
        for k in ('a', 'b'):
            best[k] = min(best[k], gpb.timed(gpb.control, repeat=1))
    control_ratio = best['a']/best['b']
    print('control ratio %.3f' % control_ratio, file=sys.stderr, flush=True)

    out = {'note': 'PREM chord at cos(theta_z) = -0.9; Earth analogue of Fig. 11',
           'costhz': ch['costhz'], 'baseline_km': ch['baseline_km'],
           'machine': platform.platform(), 'date': datetime.date.today().isoformat(),
           'control_ratio': control_ratio,
           'referee': 'prem_chord_reference.json (segment-aligned mpmath Richardson)',
           'magnus_called_via': ('osc_prob_matter_std_potential, rtol sweep, with the '
                                 'sixteen PREM layer crossings passed as t_breakpoints'),
           'cases': []}
    if OUT.exists():
        out = json.loads(OUT.read_text())

    for rc in refs['cases']:
        d = rc['flavours']
        ref = np.array(rc['reference'])
        case = next((c for c in out['cases'] if c['flavours'] == d), None)
        if case is None:
            case = {'flavours': d, 'baseline_inv_ev': ch['baseline'],
                    'energy_ev': list(refs['energy_ev']),
                    'reference': rc['reference'],
                    'reference_self_convergence': rc['self_convergence'],
                    'series': []}
            out['cases'].append(case)
        have = {s['name'] for s in case['series']}
        print('--- %d flavours ---' % d, file=sys.stderr, flush=True)

        for order in ORDERS:
            name = 'Magnus' if order == 4 else 'Magnus, order %d' % order
            if name in have:
                continue
            print('    %s' % name, file=sys.stderr, flush=True)
            pts = magnus_points(d, prof, gpb.MG_DIALS, ref, order,
                                ch['edges'][1:-1])
            case['series'].append({'name': name, 'dial': 'rtol',
                                   'magnus_exp_order': order, 'points': pts})
            OUT.write_text(json.dumps(out, indent=1))

        if 'NuOscProbExact' not in have:
            print('    NuOscProbExact', file=sys.stderr, flush=True)
            npe = gpb.npe_points(d, prof, gpb.NPE_DIALS, ref)
            if npe:
                case['series'].append({'name': 'NuOscProbExact', 'dial': 'n_slabs',
                                       'points': npe})
                for p in npe:
                    print('      n_slabs %6d  %10.2f us/prob  err %.3e'
                          % (p['n_slabs'], p['us_per_probability'],
                             p['max_abs_error']), file=sys.stderr, flush=True)
            else:
                print('      no %dnu route' % d, file=sys.stderr, flush=True)
            OUT.write_text(json.dumps(out, indent=1))

    OUT.write_text(json.dumps(out, indent=1))
    print('wrote %s' % OUT.name, file=sys.stderr)


if __name__ == '__main__':
    main()
