# Issue #92: recomputing notebook 28's cache in parallel (PR #119)

Four copies of the repository recompute disjoint sets of the value-only sections of
`notebooks/paper_figure_cache.json`; the timed `prob_vs` sections are recomputed afterwards by one
copy alone, on a quiet machine.  Only the copies are patched, never the repository.

- `nb28_par_setup.py` makes the copies and assigns the sections (`redo_only_nb28_w*.txt`,
  `redo_only_nb28_solo.txt`).
- `nb28_par_run.sh` runs the four copies, stops if any failed, waits for the machine's usual
  quiet level (1- and 5-minute load <= 1.0 for three minutes), runs the timed copy, then merges
  and diffs.  Its log is `nb28_par.log`.
- `merge_nb28.py` merges the copies' sections into one cache, in scratch.
- `diff_cache.py COMMITTED.json RECOMPUTED.json` compares two caches section by section;
  `nb28_cache_diff.txt` is its report for #92 (82 sections recomputed: 62 identical, 20 moved).

Worth rerunning when the cache has to be recomputed.  Each copy needs a `../docs` link beside it,
`MAGNUS_PAPER_REDO` reaches only `cached()`, and timings are re-measured only under
`MAGNUS_PAPER_RETIME`.
