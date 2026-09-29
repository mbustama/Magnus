# Listing 1's timing

- **`b8/`**, step B8 of #120, on `main` at 5bd8d8c.  `listing_printed.py` is the paper's Listing 1
  as printed after PR #128; `listing_default_tol.py` is the same at the default tolerance, and
  `listing_hybrid.py` the same on the hybrid strategy, for reference.  `fresh.py` times one run in
  a new process, imports apart; `warm.py` repeats it in one process and reads the engines;
  `drive.py` runs them one process at a time and writes `b8_results.json`.  At 1e-12 the four
  calls took 0.66 s from a fresh start (1.1 s with Python's start-up and the imports) and 0.48 s
  warm; at the default tolerance, 0.21 s and 0.07 s; on the hybrid, 10.4 s.
- `time_lst_validation.py`, `time_lst_default.py`, `time_lst_default_reversed.py` and their logs:
  the original fresh-process protocol, median of five.
- `run_listing.py`, `run_listing_warm.py`, `warm_both.py`: #70's versions, per curve.

To rerun B8 on a quiet machine: `python b8/drive.py PARENT`, where `PARENT/b8/` holds these files.
It overwrites `b8_results.json`, so run it on a copy.  If the listing changes, extract it again
from `resources/paper/main.tex`.
