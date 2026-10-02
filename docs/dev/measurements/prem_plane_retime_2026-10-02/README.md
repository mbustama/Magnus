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

The new series were written into `notebooks/external_earth_plane.json` (three flavors) and
`notebooks/external_prem_speed_accuracy_new.json` (`sterile_3plus1`), Magnus series only, with
the record `_magnus_provenance.retimed_series`. Notebook 28 then re-rendered the figure,
cache-only.
