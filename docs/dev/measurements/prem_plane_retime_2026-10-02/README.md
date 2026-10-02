# Fig. `prem_plane`: Magnus's Earth-plane series re-timed (2026-10-02)

Magnus's points on the two Earth speed-accuracy planes of Fig. `prem_plane`
(`resources/paper/figs/speed_accuracy_combined.pdf`) re-measured on the paper's laptop at commit
35dcdd5, after PR #210 let `notebooks/gen_prem_plane_magnus.py` run at three flavors again.
The slab-count series had been measured on 2026-09-02 at 6751a8c and was stale: Magnus is now
1.6 to 3.9 times faster on it. The tolerance series, re-timed on 2026-09-28, moved by under 11%.
The errors are unchanged except the 1-slab points and three points near the reference floor.

- `timed_run.log`: the run's printed output (one line per point; 15.5 minutes, 31 points).
  It started at load 1.28, with nothing else running, on mains power, under governor
  powersave.
- `magnus_timed.json`: the run's output, written by
  `python notebooks/gen_prem_plane_magnus.py timed /tmp/magnus_timed.json`.
- `compare_output.txt`: `python notebooks/gen_prem_plane_magnus.py compare /tmp/magnus_timed.json`,
  the old against the new time and error at every point.
- `adopt_retimed_series.py`: the code that wrote the new series into the stored files, run
  inline from the repository root and kept as it was run. Fed `magnus_timed.json` and the
  files as they were at a978203^, it reproduces both files of a978203 byte for byte.

The new series were written into `notebooks/external_earth_plane.json` (three flavors) and
`notebooks/external_prem_speed_accuracy_new.json` (`sterile_3plus1`), Magnus series only, with
the record `_magnus_provenance.retimed_series`. Notebook 28 then re-rendered the figure,
cache-only.

## Added the same day: rtol 1e-5, 1e-7 and 1e-9

The tolerance series had gaps at 1e-5, 1e-7 and 1e-9. They were timed from 16:28 to 16:33 on the
same laptop and package code (`src/` unchanged from 35dcdd5; the repository was at 0dc1e7c),
with the stored points rtol 1e-4, 1e-6 and 1e-8 interleaved as controls. The controls came
back within 1.1% of their stored times (0.989 to 1.002) with identical errors, so the new
points were written and every stored point was left as it was.

- `added_rtol_accuracy.json`: the untimed check made before the run, written by
  `python notebooks/gen_prem_plane_magnus.py accuracy OUT.json --series tolerance --knobs=-4,-5,-6,-7,-8,-9`.
  It gave the new points' errors, and the controls reproduced their stored errors exactly.
- `added_rtol_run.log`: the run's printed output (one line per point). It started at load 0.92,
  with nothing else running, on mains power, under governor powersave.
- `added_rtol_timed.json`: the run's output, written by
  `python notebooks/gen_prem_plane_magnus.py timed OUT.json --series tolerance --knobs=-4,-5,-6,-7,-8,-9`,
  with OUT.json in the session's scratch directory.
- `add_rtol_points.py`: wrote the three new points of each panel and the record
  `_magnus_provenance.added_tolerance_points`. It refuses if a control's error differs from the
  stored one or its time moved by more than 10%. Fed `added_rtol_timed.json` and the files as
  they were at 0dc1e7c, it reproduces both files of a079f56 byte for byte.

rtol 1e-9 is stored but not drawn: it reaches 4.7e-11 at three flavors, beneath that panel's
lower edge, and 2.2e-11 at 3+1, on that panel's lower edge (2e-11). The figure stops at
`MAGNUS_RTOL_FLOOR` = 1e-8 in `notebooks/make_notebooks.py`, as before.
