# Issue #192: the numbers of paper Sec. 6, re-run against the code

Each script reproduces the numbers of one stretch of Sec. 6 of `resources/paper/main.tex`
(lines as on main at 4083ad8), on one core. Timings are not checked here; they belong to #194.

| Folder | Lines | Content |
|---|---|---|
| `A/` | 1383–2422 | Parameter sets, batched calls, constant density, wrappers vs. scenario vs. direct calls, the Hamiltonians |
| `B/` | 2423–3257 | Plotting, probability vs. energy and distance, layered matter, bi-probability, the Earth, the strategy report |
| `C/` | 3258–3763 | The Sun, the supernova shock (`shockdef.py` keeps notebook 14's pre-1.2.0 normalization), astrophysical flavor |
| `D/` | 3764–4649 | Cavity, geoneutrinos, solar tomography, jet, long-range force, turbulence, `n_jobs`, the CLI |

Run any script from the repository root, e.g. `python docs/dev/measurements/issue192_sec6_numbers/B/s3_bip.py`.
Some scripts write `.npy`/`.npz` intermediates to the working directory; those are not kept.

The results (about 320 numbers confirmed, 18 wrong as written, 8 borderline) are in
https://github.com/mbustama/Magnus/issues/108#issuecomment-5934731733.
