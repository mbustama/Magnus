# Issue #122: the clamped step onto the slab cap (PR #123)

Run on `main` at 34107d9 and on the fix branch.

- `verify122.py OUT.json` records probabilities, warnings and times of 42 batched scans (Listing
  1's curves from 1e-3 to 5e-13, with slab floors; longer exponential profiles; a breakpoint grid;
  a solar scan) for one source tree.  `v122_old_1.json` (`main`) and `v122_narrow.json` (the fix
  as merged) are its outputs: identical probabilities, and the caps warning gained by the two
  scans that stopped on the clamped step.
- `clamp_cases.py` logs the agreement at the clamped step, to `clamp122.jsonl` and
  `clamp_cases_P.json`, for the #94 test's case and both #122 cases.  It ran on a copy patched
  with `clamp122_refusal_and_hook.patch`, which holds the first, blunter version of the fix and
  the hook.  The numbers are the table in `MIN_EFFECTIVE_REFINEMENT`'s docstring.
- `case94b.py` scores the #94 test's capped 108 -> 115 step against converged references: 0.49 of
  the tolerance.  That is why the merged fix keeps that step certifying.
