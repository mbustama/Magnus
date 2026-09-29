# Issues #91 and #93: re-measured series adopted into the stored benchmark files (PR #109)

Branch `adopt-71-series` from `main` at bd9f9d7, committed as 6ede520.  `RESUME.md` is that
work's handover note.  The driver behind the #93 tolerance series is committed as
`notebooks/gen_prem_plane_magnus.py`; `acc_6751.json` shows it reproducing the 16 stored errors
at 6751a8c.

- **#91:** `rerun91_prem.py` writes `prem_chord_rerun.json` and `prem_chord_rerun2.json`, the
  quiet re-measure that was adopted (control ratio 0.995; log `prem_rerun2.log`).
  `rerun91_njobs.py` writes `njobs_rerun.json`.
- **#93:** `timed93_main.json` (log `timed93_main.log`), adopted into
  `notebooks/external_earth_plane.json` and `notebooks/external_prem_speed_accuracy_new.json`;
  `acc_main_auto.json` and `acc_main_magnus.json` are the same series under both strategies.
- `adopt_91_93.py` writes the adopted series into the stored files.  `check_nb28.py`, reusable,
  reports what a notebook-28 rebuild moved against HEAD: pixels for figures, text for outputs.
- `probe93.py`, `probe93b.py`, `probe93c.py`: probes.
