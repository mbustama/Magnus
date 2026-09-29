**DONE 2026-09-28: committed as 6ede520, PR #109.**

# Resume: adopt the re-measured #91/#93 series and re-render (paused 2026-09-28)

## State
- Branch `adopt-71-series` (from `main` at `bd9f9d7`), NOT committed. Do not switch branches before resuming,
  or ask the author for a local WIP commit first.
- Done and in the working tree:
  - `notebooks/external_prem_chord_benchmarks.json`: 2nu and 3nu Magnus series (orders 4/6/8) replaced by the
    quiet re-measure (`prem_chord_rerun2.json`, control ratio 0.995); record `magnus_rerun_71` added.
  - `notebooks/external_earth_plane.json` (3nu) and `external_prem_speed_accuracy_new.json` (3+1 panel):
    the `Magnus (tolerance)` points replaced by `timed93_main.json`; `_magnus_provenance.retimed_tolerance_series`
    added. Slab-count series untouched. The 1e-3 point of each is flagged `best_at_this_accuracy = False`
    (same error as 1e-2, slower: emit_figures.py's rule), so the figure skips it.
  - `notebooks/gen_prem_plane_magnus.py` (staged): the reconstructed driver; at 6751a8c it reproduces all 16
    stored tolerance-series errors (`acc_6751.json`); on main `magnus` and `auto` give the same engine.
  - Registered in `tests/test_file_tree.py` and `docs/source/installation.rst` (file-tree test passes).
- The notebook-28 rebuild was STOPPED midway and its outputs restored to `main`; re-run it.

## Remaining (about 10 minutes)
1. `MAGNUS_PAPER_CACHE_ONLY=1 python notebooks/make_notebooks.py --only 28` (~2 min).
2. `python resources/handover_auto_tight/adopt_71/check_nb28.py`. Expected to change: `smooth_reach.pdf` (#91)
   and `speed_accuracy_combined.pdf` (#93). Restore every pixel-identical PDF from HEAD, and
   `phase_vs_profile.pdf` too (its top-right panel is live round-off noise). Look at the two changed figures.
   Keep the notebook; keep gallery PNGs only if they differ in pixels.
3. Gates: `pytest tests/test_notebooks_match_their_generator.py tests/test_file_tree.py $(grep -l "28_magnus\|img/gallery\|paper/figs" tests/*.py) -q`;
   `ruff check src/magnus/ tests/ docs/ notebooks/make_notebooks.py`;
   strict docs into a fresh dir (installation.rst changed): `sphinx-build -n -W --keep-going -b html docs/source <scratch>/docs_strict_N`.
4. STOP: report, then commit/push/PR only on the author's word. The PR closes #91 and #93. The paper edits are
   in #108 (items 5-6 hold after the re-measure; items 13-14 come from #93). Overleaf then needs
   `smooth_reach.pdf` and `speed_accuracy_combined.pdf`.
5. Still open afterwards: #92 (full notebook-28 recompute in a scratch copy, hours, overnight); Part B of
   `../HANDOVER.md` (auto takes the ladder at tight tolerances).

## Evidence in this folder
`timed93_main.json/.log` (the #93 timed pass: 9 min, worst block cv 0.017), `prem_chord_rerun*.json`,
`prem_rerun2.log`, `njobs_rerun.json` (n_jobs=1 shared arm: 0.111 -> 0.114 s, no text change), `acc_*.json`
(untimed accuracy at 6751a8c and on main), `adopt_91_93.py` (the adoption), `probe93*.py`, `rerun91_*.py`.
