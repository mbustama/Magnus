# -*- coding: utf-8 -*-
"""Acceptance test for the cross-method agreement diagnostic (robustness programme, item 1).

The question this answers is NOT "do the engines agree on today's code" -- they do, and that
proves nothing about the instrument.  It is:

    **on the constructions where a method was silently wrong BEFORE the fixes, would a
    cross-check between engines have seen it?**

So this script is written to run against *either* tree, and is driven by ``PYTHONPATH``:

```bash
git worktree add /tmp/mainwt 978663a
PYTHONPATH=/tmp/mainwt/src python crosscheck_acceptance.py > cc_prefix.txt
python crosscheck_acceptance.py > cc_fixed.txt
```

It deliberately does **not** import ``oscprob.cross_check_strategies``: that function does not
exist on the pre-fix tree, and an acceptance test that can only run where the fix already is
would be testing nothing.  The engine forcing here is the battery-3 spy pattern -- monkeypatch
the dispatchers, read back which one answered -- and reproduces what
``cross_check_strategies`` does through ``_ENGINES_DISABLED``.

PASS CRITERION, stated before running, in a form that does not depend on which tree it runs
against -- a criterion phrased as "these seven must be detected" would pass on the pre-fix tree
and fail on the fixed one for the best possible reason (the defects are gone), which is not a
useful instrument:

  (X1) **Whenever at least one engine is outside the requested tolerance, the maximum
       cross-family spread must be at least that tolerance.**  A cross-check detects
       disagreement, so it must see a wrong answer whenever some *other* engine got it right.

  (X2) The one exception is stated in advance rather than discovered afterwards: when
       **every** engine is outside tolerance, they are wrong together and there is nothing to
       disagree about.  ``FINDINGS §8.3`` records exactly one such construction -- a Gaussian
       narrower than the probe spacing, invisible to every grid the package lays down -- and it
       is carried here as a row that is *expected* to go undetected.  Reporting that honestly
       is the point; a diagnostic that claimed to cover it would be lying.

On the pre-fix tree X1 exercises all seven ``FINDINGS §8.2`` silent misses, which is the
acceptance the robustness brief asks for.  On the fixed tree most of those rows have no wrong
engine left, so X1 is vacuous for them and the run doubles as a regression check.
"""

import sys
import warnings

import numpy as np

import harness as H
import magnus.oscprob as oscprob
from battery2 import bump_profile, many_bumps_profile, ne_res_for

TOL = 1e-3
L0, L1 = 0.0, 1.0*H.L_SCALE
SPAN = L1 - L0


# ---------------------------------------------------------------- engine forcing
# (strategy, cumulative-or-None, engines to switch off so that this one is reached)
FORCING = {
    'hybrid':     ('hybrid', None, ('ip_exp', 'separable')),
    'ip_exp':     ('magnus', None, ('hybrid', 'separable')),
    'separable':  ('magnus', None, ('hybrid', 'ip_exp')),
    'magnus':     ('magnus', None, ('hybrid', 'ip_exp', 'separable')),
    'cumulative': ('magnus', True, ('hybrid', 'ip_exp', 'separable')),
}

FAMILY = {'hybrid': 'adiabatic', 'ip_exp': 'interaction-picture',
          'magnus': 'magnus-ladder', 'cumulative': 'magnus-ladder',
          'separable': 'magnus-ladder'}

_DISPATCHERS = {'hybrid': '_osc_prob_hybrid_dispatch',
                'ip_exp': '_osc_prob_ip_exp_dispatch',
                'separable': '_osc_prob_scan_separable_dispatch'}


class Forced:
    """Switch off the named dispatchers and record which engine answered.

    Same shape as ``battery3.Spy``, extended to disable as well as observe.  Restores every
    patched name on the way out, including on an exception -- ``cumulative=True`` raises by
    design on a request it cannot serve.
    """

    def __init__(self, off=()):
        self.off = tuple(off)
        self.answered = set()

    def __enter__(self):
        self._saved = {}
        for label, name in _DISPATCHERS.items():
            orig = getattr(oscprob, name)
            self._saved[name] = orig

            def make(label=label, orig=orig):
                def wrapper(*a, **k):
                    if label in self.off:
                        return NotImplemented
                    r = orig(*a, **k)
                    if r is not NotImplemented:
                        self.answered.add(label)
                    return r
                return wrapper

            setattr(oscprob, name, make())
        cum = oscprob._osc_prob_cumulative_scan
        self._saved['_osc_prob_cumulative_scan'] = cum

        def cum_wrapper(*a, **k):
            self.answered.add('cumulative')
            return cum(*a, **k)

        oscprob._osc_prob_cumulative_scan = cum_wrapper
        return self

    def __exit__(self, *e):
        for name, orig in self._saved.items():
            setattr(oscprob, name, orig)
        return False


def run_engine(label, ne, d, params, energy, baseline):
    """One engine's answer, or (None, reason)."""
    strategy, cumulative, off = FORCING[label]
    kw = dict(strategy=strategy)
    if cumulative is not None:
        kw['cumulative'] = cumulative
    try:
        with Forced(off) as f, warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            P = np.asarray(oscprob.osc_prob_matter_std_potential(
                d, ne, energy, baseline, params, L0=L0,
                density_is_of_number_of_electrons=True, **kw))
    except Exception as exc:                       # noqa: BLE001 -- a decline, not a failure
        return None, type(exc).__name__
    # 'magnus' is the terminal path: there is no dispatcher to spy on, so it answered exactly
    # when nothing else did.
    answered = f.answered if label != 'magnus' else (f.answered or {'magnus'})
    if label not in answered:
        return None, 'declined (answered by %s)' % (','.join(sorted(f.answered)) or 'magnus')
    names = sorted({w.category.__name__ for w in caught})
    return (P.reshape(d, d), ','.join(names) or '-')


def cross_check(ne, d, params, energy, baseline):
    answers, notes = {}, {}
    for label in FORCING:
        P, note = run_engine(label, ne, d, params, energy, baseline)
        if P is None:
            notes[label] = note
        else:
            answers[label], notes[label] = P, note
    spread, best, best_pair = {}, 0.0, None
    labs = list(answers)
    for i, a in enumerate(labs):
        for b in labs[i + 1:]:
            s = float(np.max(np.abs(answers[a] - answers[b])))
            spread[(a, b)] = s
            if FAMILY[a] != FAMILY[b] and s > best:
                best, best_pair = s, (a, b)
    return answers, notes, spread, best, best_pair


# ---------------------------------------------------------------- constructions
def build_cases():
    """The FINDINGS §3 constructions, verbatim from the batteries that produced them."""
    p2 = H.params_for(2)
    lo, hi = 0.02*H.NE0, 0.30*H.NE0
    mid = 0.5*L1
    ner2 = ne_res_for(2, p2, 10.0e6)
    rng = np.random.default_rng(7)
    lc = L0 + (0.37 + 0.2*rng.random())*SPAN         # battery2's pinned "random position"

    def step(l):
        x = np.asarray(l, dtype=float)
        return H.scalarize(np.where(x < mid, lo, hi))

    def kink(l):
        x = np.asarray(l, dtype=float)
        return H.scalarize(lo + (hi - lo)*np.abs(x - mid)/L1)

    def near_sing(l):
        x = np.asarray(l, dtype=float)
        return H.scalarize(lo + (hi - lo)*0.001/(np.abs(x - 1.02*L1)/L1 + 1e-3))

    def sine_span_over_7(l):
        per = SPAN/7.0
        return H.scalarize(ner2*(1.0 + 0.9*np.sin(2.0*np.pi*np.asarray(l, float)/per)))

    ten = many_bumps_profile(ner2, L0 + SPAN*(np.arange(10) + 0.5)/10, 0.25*SPAN/10)

    # (label, ne, d, params, energy, baseline, prefix-tree error from FINDINGS §8.2/§3,
    #  whether it was SILENT there)
    return [
        ('step function, unmarked edge', step, 2, p2, 50.0e6, L1, 5.395e-01, True),
        ('kink, C0 but not C1', kink, 2, p2, 50.0e6, L1, 1.448e-02, True),
        ('singularity approached', near_sing, 2, p2, 50.0e6, L1, 8.625e-03, True),
        ('sinusoid, period span/7', sine_span_over_7, 2, p2, 10.0e6, L1, 1.672e-02, True),
        ('ten crossings', ten, 2, p2, 10.0e6, L1, 3.907e-02, True),
        ('sub-threshold bump, w=3e-2 span', bump_profile(ner2, lc, 3e-2*SPAN), 2, p2,
         10.0e6, L1, 4.388e-03, True),
        ('sub-threshold bump, w=1e-2 span', bump_profile(ner2, lc, 1e-2*SPAN), 2, p2,
         10.0e6, L1, 7.701e-03, True),
        ('narrow bump, w=3e-5 span (§8.3)', bump_profile(ner2, lc, 3e-5*SPAN), 2, p2,
         10.0e6, L1, 2.907e-02, False),
    ]


def main():
    print('# Cross-check acceptance.  magnus from: %s' % oscprob.__file__)
    print('# requested tolerance %.0e; L0=%.3e L1=%.3e\n' % (TOL, L0, L1))
    verdicts = []
    for (label, ne, d, params, energy, baseline, prefix_err, was_silent) in build_cases():
        H_of_l = H.H_factory(d, params, H.vcc_of(ne), energy)
        Pref = H.P_of(H.exact_U(H_of_l, L0, baseline, d))
        answers, notes, spread, best, best_pair = cross_check(ne, d, params, energy, baseline)

        print('%s   [FINDINGS: %.3e%s]' % (label, prefix_err,
                                           ', SILENT' if was_silent else ', all engines wrong'))
        for lab in FORCING:
            if lab in answers:
                print('    %-11s err=%9.3e   warns=%s'
                      % (lab, H.maxabs(answers[lab] - Pref), notes[lab]))
            else:
                print('    %-11s %s' % (lab, notes[lab]))
        for (a, b), s in sorted(spread.items(), key=lambda kv: -kv[1]):
            print('    spread %-11s vs %-11s %9.3e%s'
                  % (a, b, s, '' if FAMILY[a] != FAMILY[b] else '   (same family)'))

        errs = {lab: H.maxabs(P - Pref) for lab, P in answers.items()}
        any_wrong = max(errs.values()) > TOL
        all_wrong = min(errs.values()) > TOL
        detected = best >= TOL
        if all_wrong:
            state, ok = 'every engine wrong -- X2, not detectable', True
        elif any_wrong:
            state, ok = 'an engine is wrong -- X1 requires detection', detected
        else:
            state, ok = 'no engine outside tolerance -- nothing to detect', True
        print('    -> max cross-family spread %.3e via %s\n       %s   %s\n'
              % (best, best_pair, state,
                 'ok' if ok else '<-- X1 FAILURE: a wrong engine went unseen'))
        verdicts.append((label, any_wrong, all_wrong, best, detected, ok))

    print('=== SUMMARY ===')
    testable = [v for v in verdicts if v[1] and not v[2]]
    print('rows where some engine is outside %.0e : %d   (X1 applies)' % (TOL, len(testable)))
    print('    of which detected                  : %d' % len([v for v in testable if v[4]]))
    print('rows where EVERY engine is wrong        : %d   (X2, expected undetectable)'
          % len([v for v in verdicts if v[2]]))
    print('X1 failures                             : %d'
          % len([v for v in verdicts if not v[5]]))
    for label, any_w, all_w, best, det, ok in verdicts:
        tag = ('every engine wrong' if all_w
               else ('detected' if any_w and det
                     else ('MISSED' if any_w else 'all engines inside tolerance')))
        print('   %-34s spread %9.3e  %s' % (label, best, tag))


if __name__ == '__main__':
    sys.exit(main())
