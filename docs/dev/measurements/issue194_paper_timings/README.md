# Issue #194: the paper's timings, re-measured on one quiet machine

`run_all.py` reads each timed listing or snippet out of `resources/paper/main.tex` (by its
`\label` or a unique anchor string), runs it as printed, and times it: warm-up excluded, median
of N repeats with min/max, fresh processes where the paper quotes one. Items where the paper
gives a number but prints no code (2, 3, 5b, the eight arrangement curves) are built from the
text, and their rows say so. Threads are pinned to one per process (`OMP_`, `MKL_`, `OPENBLAS_`,
`NUMEXPR_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`), as in `notebooks/gen_njobs_scaling.py`.

To run, from the repository root, on the laptop with nothing else running (plugged in,
no browser): first `python docs/dev/measurements/issue194_paper_timings/run_all.py --quick`
(about 2-3 minutes; checks the script works, the numbers mean nothing), then the real run:

    python docs/dev/measurements/issue194_paper_timings/run_all.py

This takes about 5 minutes. `--only 1,4,7` reruns some of the items (1-9, as in #194).

Send back the whole printed output, plus the file `results_<hostname>_<date>.json` written
next to the script. The table gives the line in `main.tex`, what the paper says, the
measurement, and suggested text rounded as the paper rounds. Check the stability-control
line first: if the drift ratio goes above 1.10, the machine was busy and that run should be
repeated. The last lines print the command that re-times notebook 28's cache (#194 part 2).
