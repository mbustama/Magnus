# Issue #120: `strategy='auto'` takes the ladder at tight tolerances (PR #128)

The measurement (B1) is `tools/auto_tight/`.  Its output on `main` at 34107d9 is
`tight_b1.jsonl`, from which `AUTO_LADDER_TIGHT_MAX_PHASE`'s docstring numbers come.  The rest is
the double check of the rule before it was implemented, on a copy of `main` at 6b24be1 patched
with `route120_rule_and_hook.patch`: the rule, plus a log of every decision beside the one `main`
makes.

- `b2_table.py tight_b1.jsonl` applies the rule to B1's records under each reading of the
  tolerance, with and without the 2 000 rad cap.  `cmp_b1.py` compares B1 before and after #122
  (identical probabilities, warnings and engines).  `b1_margins.py` measures the /2 and /3
  margins.
- `route120_tests.jsonl` and `route120_nb.jsonl` are every routing decision in the test suite and
  in the 29 notebooks on the patched copy, and `route120_report.py` summarizes them: none of 1238
  changes at 1e-6 and looser, and below it 5 change in the suite and none in the notebooks.
  The final code replays all of them with no mismatch.
- `retime120.py` re-times, alone and interleaved, the four routed cases once slower than the
  hybrid (0.38 to 0.54 of its time).
- `scans120.py`: two energy scans outside B1 (25 and 250 km).  `table_rows.py`: the rows of
  `engines.rst`'s engine table, on either tree.  `probe_tests.py`, `probe_tests2.py`: what the
  tests of PR #128 assert, measured before they were written.

Listing 1's timing for the paper is in `../listing1_timing/b8/`.
