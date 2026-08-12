# Interaction sweeps

Five scripts kept because the defects they found were **interactions between features**,
not defects in any one feature, and nothing else in the repository looks for that shape.

This package's history is the argument. A 7440x accuracy hole survived because each axis
was swept alone — clustered spectra at norm 1, and large norms at generic separation — and
the damage needed both at once. The `A2b` sterile matter term was invisible at three
flavours. The Earth's layered `Y_e` and the sterile projector each behaved correctly and
described different media when combined.

Run them from the repository root. None needs an idle machine: every one is an accuracy or
structural check, not a timing measurement.

| script | what it crosses | what it found |
|---|---|---|
| `check_angles.py` | signature vs docstring vs body, for `angles` | four functions that accepted `angles` and silently ignored it — invisible to tests, docs build and ruff |
| `check_defaults.py` | every documented `Default: X` against the signature | nothing outstanding; one known false positive (`2*pi` vs `6.28318…`) |
| `cross_features.py` | flavours x environment x `nubar` x `average` x the four `angles` conventions, 36 cells | clean: worst convention spread 8.3e-13, worst unitarity 1.7e-12 |
| `cross_features2.py` | `strategy`, `expm_backend`, the antineutrino identity, `integration_method` | clean; taught that `strategy` is correctly refused in vacuum |
| `cross_features3.py` | `cumulative`, `n_jobs`, `t_breakpoints`, the per-layer `electron_fraction_*` overrides against the O50 guard | `t_breakpoints` raised `TypeError` on every Earth wrapper; `n_jobs` moves the answer at the requested tolerance |

Two of these have since been promoted into the suite and run in CI:
`tests/test_angles.py` carries the structural checks and the guards, and
`tests/test_routine_listings.py` the module-listing check. The sweeps here are the broader,
slower exploration those were distilled from — keep them for the next audit rather than
re-deriving the grid.

**A caution learned twice while writing them.** When a sweep disagrees, suspect the harness
first. `cross_features2.py` initially demanded that `strategy` work in vacuum, where there
is no matter potential and the library refuses it with an explicit message; `cross_features3.py`
initially demanded bitwise `n_jobs` agreement, which was never the contract. Both times the
library was right. The comments in those files record the corrected expectation rather than
hiding the mistake, because the wrong expectation is the more instructive half.
