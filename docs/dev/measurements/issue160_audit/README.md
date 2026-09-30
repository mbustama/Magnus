# Issue #160: every sub-case of the checklist, run

`audit160.py` runs one case per sub-case of issue #160's checklist -- 733 of them -- and reports
each that does not behave as the checklist asks: a refusal that names the argument, a warning,
a documented rule, or a call that must still work.  Every group opens with a `BASE` case that must
pass, so a harness call that is itself broken cannot pass by raising.

It backs the "Closes" line of the bundled pull request for #155, #160 and the issues folded into
it, and the reproducer table in that pull request's description.

Run it from the repository root against a source tree:

    python docs/dev/measurements/issue160_audit/audit160.py src out.json

| Tree | Cases | Pass | Fail |
|---|---|---|---|
| `main` at `0939598` (after PR #173) | 733 | 629 | 104 |
| the bundled branch | 733 | 733 | 0 |

Nine cases are marked "Revised" in the script: they expect the behaviour decided on the branch
rather than the one the checklist first asked for, and each says why.

* `cumulative=True` on a constant density is not ignored: it runs the cumulative engine.
* `rtol` in vacuum stays accepted, since the command line and shared calls forward one set of
  numerics to every wrapper; the vacuum docstrings say it has no effect there.
* A points-per-slab setting under `integration_method='gl'` warns
  (`IgnoredQuadratureSettingWarning`) and is documented, rather than refused.
* `adiabatic.oscillation_sampling` and `find_hidden_features` return a quiet empty report for a
  broken profile rather than raising: their documented design, pinned by
  `tests/test_adiabatic.py::test_oscillation_sampling_refuses_quietly_rather_than_breaking_the_call`.

The per-fix tests, each failing on `main` and passing on the branch, are in
`tests/test_validation_160_reopened.py`.
