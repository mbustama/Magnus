r"""Magnus's series on the two Earth speed-accuracy planes of Fig. ``prem_plane`` (issue #93).

The Magnus series in ``external_earth_plane.json`` (three flavors) and in the 3+1 panel of
``external_prem_speed_accuracy_new.json`` were first measured on 2026-09-02 (commit 6751a8c) by a
driver kept outside the repository.  This is its reconstruction: at 6751a8c it reproduces the
errors that run stored.  It registers Magnus with the benchmark harness at run time, so
``resources/benchmarks`` stays byte-identical to its source.

The problem, as ``magnus_own_reference.json`` defines it:

* the PREM chord at cos(theta_z) = -0.9, electron fraction 0.5 and neutron-to-proton ratio 1,
  with twelve energies: 3 to 40 GeV at three flavors (the manifest's CHORD/12x1 grid), 300 GeV
  to 30 TeV at 3+1;
* the manifest's oscillation parameters, through the harness's ``conversions.py``, plus
  sin^2 theta_14 = sin^2 theta_24 = 0.1, theta_34 = 0 and Delta m^2_41 = 1 eV^2 at 3+1;
* the knob: ``n_slabs = knob`` with no tolerance (knob > 0), or ``rtol = 10**knob`` with
  ``atol = rtol/100`` (knob < 0);
* one call per evaluation, for the whole probability matrix at the twelve energies, scored as
  the nu_mu row at three flavors and as P(nu_mu -> nu_mu) at 3+1: the largest absolute deviation
  from Magnus's own 50-digit reference.

Timing is the harness's AMORTIZED protocol, ``runner.amortized``: a 25-step delta_CP scan,
autoranged to blocks of at least 0.25 s, 30 to 100 blocks.  ``us_per_probability`` is the mean
over the blocks, per energy, as the stored series record it.

From the repository root:

    python notebooks/gen_prem_plane_magnus.py accuracy OUT.json [--series ...] [--strategy ...]
    python notebooks/gen_prem_plane_magnus.py timed OUT.json [--series ...] [--strategy ...]
    python notebooks/gen_prem_plane_magnus.py compare OUT.json

``accuracy`` is untimed; ``timed`` adds the harness's timing to every point; ``compare`` prints
OUT.json against the stored series.  Only OUT.json is written: the stored files change on the
author's decision alone.  No clock is named in this file, since the harness refuses an adapter
whose source names one.
"""
import argparse
import datetime
import json
import math
import os
import pathlib
import platform
import subprocess
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'resources'/'benchmarks'/'bench'))
sys.path.insert(0, str(HERE))

import runner                                                      # noqa: E402
import magnus.earth as earth                                       # noqa: E402
import magnus.globaldefs as gd                                     # noqa: E402
import magnus.oscprob as oscprob                                   # noqa: E402

REFERENCE = json.loads((HERE/'magnus_own_reference.json').read_text())
COSTHZ = -0.9
PANELS = {'3nu': dict(reference='three_flavor', n_flavors=3,
                      stored=('external_earth_plane.json', None),
                      slabs=(1, 2, 4, 8, 16, 32, 64, 128, 256)),
          '3+1': dict(reference='sterile_3plus1', n_flavors=4,
                      stored=('external_prem_speed_accuracy_new.json', 'sterile_3plus1'),
                      slabs=(1, 8, 32, 128, 512, 2048))}
TOLERANCE_KNOBS = (-1, -2, -3, -4, -6, -8, -10, -12)
AMORTIZED = dict(samples=30, steps=25, min_block=0.25, max_samples=100)
CHANNELS = {3: ('numu->nue', 'numu->numu', 'numu->nutau'), 4: ('numu->numu',)}
STRATEGY = 'magnus'


class MagnusAdapter(object):
    r"""Magnus as one of the harness's codes; see the module docstring."""
    name = 'Magnus'

    def capabilities(self):
        return {'batches_energy': True, 'batches_zenith': False,
                'batch_symbol': ('oscprob.osc_prob_3nu_earth / osc_prob_4nu_earth over the '
                                 'energy array, nu_i = nu_f = None'),
                'knob_name': 'n_slabs=knob (knob>0) | rtol=10^knob, atol=rtol/100 (knob<0)',
                'channels': list(CHANNELS[self._n])}

    def environment(self):
        return {'magnus': oscprob.__file__, 'strategy': self._strategy}

    def setup(self, problem):
        self._n = getattr(problem, 'n_flavors', 3)
        # From the problem: the harness imports a second copy of this module, whose STRATEGY
        # never sees the command line.
        self._strategy = getattr(problem, 'strategy', STRATEGY)
        self._e_ev = np.asarray(problem.energies_gev, dtype=float)*gd.UNIT_GEV
        par = REFERENCE['oscillation_parameters']
        self._fixed = dict(
            costhz=COSTHZ, L=earth.distance_traveled_inside_earth(COSTHZ)*gd.CONV_KM_TO_INV_EV,
            s12=math.sqrt(problem.s12sq), s23=math.sqrt(problem.s23sq),
            s13=math.sqrt(problem.s13sq), D21=problem.dm21, D31=problem.dm31,
            electron_fraction=problem.ye, ratio_number_neutrons_to_protons=1.0,
            nu_i=None, nu_f=None, strategy=self._strategy)
        if self._n == 4:
            self._fixed.update(s14=math.sqrt(par['sinsq_th14']), s24=math.sqrt(par['sinsq_th24']),
                               s34=math.sin(par['th34']), D41=par['dmsq41_ev2'])
        knob = int(problem.knob)
        self._fixed.update(dict(n_slabs=knob, rtol=None, atol=None) if knob > 0
                           else dict(rtol=10.0**knob, atol=10.0**knob*1.0e-2))
        self._call = oscprob.osc_prob_3nu_earth if self._n == 3 else oscprob.osc_prob_4nu_earth
        self.configure(problem.dcp)

    def configure(self, v):
        self._dcp = v

    def reset(self):
        r"""Nothing to reset: every call rebuilds the Hamiltonians and the slab product."""

    def _row(self):
        m = np.asarray(self._call(self._e_ev, dCP=self._dcp, **self._fixed), dtype=float)
        cols = (gd.NUE, gd.NUMU, gd.NUTAU) if self._n == 3 else (gd.NUMU,)
        return m[:, gd.NUMU, :][:, cols]                  # (energy, channel)

    def evaluate(self):
        return float(np.sum(self._row()[:, -1 if self._n == 4 else 1]))

    def probabilities(self):
        return [float(v) for v in self._row().reshape(-1)]


ADAPTER = MagnusAdapter


def problem(panel, knob):
    spec = PANELS[panel]
    ref = REFERENCE[spec['reference']]
    p = runner.build_problem('CHORD/12x1', knob, 0, 0)
    p.energies_gev = list(ref['energy_gev'])
    p.n_flavors = spec['n_flavors']
    p.ye = REFERENCE['electron_fraction']
    p.strategy = STRATEGY
    par = REFERENCE['oscillation_parameters']
    got = (p.s12sq, p.s13sq, p.s23sq, math.degrees(p.dcp), p.dm21, p.dm31)
    want = (par['s12sq'], par['s13sq'], par['s23sq'], par['dcp_deg'], par['dmsq21_ev2'],
            par['dmsq31_ev2'])
    if not np.allclose(got, want, rtol=1e-12, atol=0.0):
        raise SystemExit('the manifest parameters %s differ from the reference\'s %s' % (got, want))
    return p


def driver_for(p):
    runner.CODES['magnus'] = pathlib.Path(__file__).stem
    driver = runner.load_adapter('magnus')        # refuses an adapter that names a clock
    driver.setup(p)
    return driver


def score(panel, driver):
    spec = PANELS[panel]
    ref = REFERENCE[spec['reference']]
    want = np.asarray(ref['numu_row'] if spec['n_flavors'] == 3 else ref['p_numu_numu'],
                      dtype=float).reshape(len(ref['energy_gev']), -1)
    got = np.asarray(driver.probabilities(), dtype=float).reshape(want.shape)
    return float(np.max(np.abs(got - want)))


def point(panel, knob, timed):
    p = problem(panel, knob)
    driver = driver_for(p)
    driver.configure(p.scan_base())
    info = {}
    driver._fixed['strategy_info'] = info
    err = score(panel, driver)
    del driver._fixed['strategy_info']
    out = dict(label=('%.0e' % 10.0**knob) if knob < 0 else str(knob), knob=knob,
               rtol=10.0**knob if knob < 0 else None, n_slabs=knob if knob > 0 else None,
               max_abs_error=err, engine=info.get('engine'))
    if timed:
        stats, _ = runner.amortized(driver, p, AMORTIZED['samples'], AMORTIZED['steps'],
                                    AMORTIZED['min_block'], AMORTIZED['max_samples'])
        out.update(us_per_probability=stats['mean'], us_sd=stats['sd'], us_min=stats['min'],
                   block_cv=stats['sd']/stats['mean'], blocks=stats['n'],
                   block_reps=stats['reps'])
    return out


def stored_series(panel):
    name, key = PANELS[panel]['stored']
    data = json.loads((HERE/name).read_text())
    series = (data[key] if key else data)['series']
    return {s['name']: s['points'] for s in series if s.get('name', '').startswith('Magnus')}


def run(mode, out_path, which, strategy):
    global STRATEGY
    STRATEGY = strategy
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=HERE,
                            capture_output=True, text=True).stdout.strip()
    record = dict(generated_by='notebooks/gen_prem_plane_magnus.py', mode=mode,
                  magnus=oscprob.__file__, repository_head=commit, strategy=STRATEGY,
                  date=datetime.date.today().isoformat(), machine=platform.platform(),
                  load_average_at_start=os.getloadavg(), manifest_sha256=runner.manifest_sha(),
                  thread_environment=runner.thread_environment(),
                  protocol=dict(AMORTIZED, name='AMORTIZED') if mode == 'timed' else 'ACCURACY',
                  panels={})
    for panel in PANELS:
        series = {}
        if which in ('tolerance', 'both'):
            series['Magnus (tolerance)'] = TOLERANCE_KNOBS
        if which in ('slabs', 'both'):
            series['Magnus'] = PANELS[panel]['slabs']
        record['panels'][panel] = {}
        for name, knobs in series.items():
            pts = []
            for knob in knobs:
                pts.append(point(panel, knob, mode == 'timed'))
                q = pts[-1]
                print('%-4s %-19s %6s  err %.3e%s  %s' % (
                    panel, name, q['label'], q['max_abs_error'],
                    ('  %9.1f us/prob (cv %.3f, %d blocks)' % (q['us_per_probability'],
                     q['block_cv'], q['blocks'])) if 'us_per_probability' in q else '',
                    q['engine']), flush=True)
                record['panels'][panel][name] = pts
                pathlib.Path(out_path).write_text(json.dumps(record, indent=1))
    print('wrote %s' % out_path)


def compare(path):
    new = json.loads(pathlib.Path(path).read_text())
    print('%s at %s (%s), strategy %s' % (new['mode'], new['repository_head'], new['date'],
                                          new['strategy']))
    for panel, series in new['panels'].items():
        old = stored_series(panel)
        for name, pts in series.items():
            by_knob = {q['knob']: q for q in old.get(name, [])}
            for q in pts:
                s = by_knob.get(q['knob'], {})
                line = '%-4s %-19s %6s  err %.3e -> %.3e' % (
                    panel, name, q['label'], s.get('max_abs_error', float('nan')),
                    q['max_abs_error'])
                if 'us_per_probability' in q and 'us_per_probability' in s:
                    line += '   us/prob %8.1f -> %8.1f (x%.2f)' % (
                        s['us_per_probability'], q['us_per_probability'],
                        q['us_per_probability']/s['us_per_probability'])
                print(line)


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('mode', choices=('accuracy', 'timed', 'compare'))
    ap.add_argument('out')
    ap.add_argument('--series', default='both', choices=('tolerance', 'slabs', 'both'))
    ap.add_argument('--strategy', default=STRATEGY)
    args = ap.parse_args()
    if args.mode == 'compare':
        compare(args.out)
    else:
        run(args.mode, args.out, args.series, args.strategy)


if __name__ == '__main__':
    main()
