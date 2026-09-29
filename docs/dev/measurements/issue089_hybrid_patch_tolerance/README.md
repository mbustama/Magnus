# PR #89: the hybrid strategy's patches converge to a tenth of a tight tolerance

Part A of `docs/dev/HANDOVER_AUTO_TIGHT.md`, step 4, run on the fix branch and on an export of
`main` before it.

- `verify_a.py` writes `a_head.npz` and `a_fix.npz`: every hybrid-path result at
  atol + rtol >= 1e-6 on both trees, for a bitwise comparison.
- Step 4b's case list, callables and DOP853 reference code are `tools/auto_tight/verify_b.py`,
  and its references on the fix are `tools/auto_tight/b_fix.jsonl`.  `b_head.jsonl` here is the
  same run on `main`, before the fix (logs `b_fix.log`, `b_head.log`).
- The inventory of the calls that reach the hybrid at a tight tolerance: `hp_hook.py` wraps
  `adiabatic.hybrid_propagator`, `hp_plugin.py` does the same as a pytest plugin, and
  `run_nb.py` with `ipy/` covers the notebooks.  Results: `hp_tests.log`, `hp_notebooks.log`;
  the hook's load test: `hp_hooktest.log`.
