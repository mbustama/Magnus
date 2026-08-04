# Adversarial-validation batteries

The scripts behind `../FINDINGS_ADVERSARIAL_VALIDATION.md`. They are kept so that every
failure in that document is **reproducible**, which the handover argues is worth more than a
silent patch. They are diagnostics, not tests: they are not collected by pytest, they print
tables rather than asserting, and several take tens of minutes.

Run them from this directory (they import each other by module name):

```bash
python battery2.py 1
```

| file | what it does |
|---|---|
| `harness.py` | oracles (`solve_ivp`/DOP853, `expm`), density profiles, warning capture |
| `harness6.py` | std / NSI / LIV Hamiltonian builders at d = 2…5, mirroring the wrappers |
| `battery2.py` | the detector's fixed probe grid (arg: sub-test `1`,`2`,`3`,`4`,`5`,`7`) |
| `battery3.py` | routing seams (`1`…`5`) |
| `battery4.py` | extreme numerics (`1`…`5`) |
| `battery5.py` | designed to break; `4 <n>` runs the n-case fuzzer |
| `battery6.py` | flavor count as a first-class axis (`1`…`7`) |
| `battery7.py` | cross-module and oracle diversity |
| `verify_b2.py` | proves the Battery 2 failures are code bugs, against four oracles |
| `diag_gamma.py` | dense γ scan separating *detection miss* from *sub-threshold accumulation* |
| `user_impact.py` | whether the Battery 2 defects reach a user through the public API |
| `bitident.py` | dumps the bit-identity set; run under both trees, then diff the `.npz` |
| `attribute.py` | same failing cases under `main` and the branch, to attribute each defect |

Added by the robustness programme (`../HANDOVER_ROBUSTNESS_PROGRAMME.md`):

| file | what it does |
|---|---|
| `crosscheck_acceptance.py` | acceptance test for `oscprob.cross_check_strategies`: would a cross-check between engines have caught the known silent misses? Runs against **either** tree |
| `invariants.py` | the oracle-free invariants (item 4), swept over a profile matrix; sets the bounds in `tests/test_invariants.py` |
| `warn_fp.py` | every warning's true- and false-positive rate, plus whether `MagnusConvergenceWarning` fires on the level whose answer is returned |
| `resolution_fp.py` | the resolution test's false-positive rate swept over **sub-intervals**, the axis the original `RESOLUTION_RATIO` measurement did not have |
| `constants_audit.py` | provenance for `fd_step_frac` (against the analytic `dH/dl`) and `threshold0` (cost against accuracy, at three tolerances) |
| `constants_audit2.py` | the seven remaining constants, over 18 workloads spanning points, baseline scans **and** energy scans |
| `weak_band.py` | where the hybrid path's self-certification is weak: window count and gamma margin against the error actually made |
| `crosscheck_benefit.py` | whether a default-path cross-check would earn its cost. It does not — see `FINDINGS_ROBUSTNESS_PROGRAMME.md` §11.2 |
| `fallback_quality.py` | every applicable engine scored on the same request. This is what moved `HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` from 25 to 8 |

Pre-existing scripts the original table omitted, listed here so the directory is fully covered:

| file | what it does |
|---|---|
| `battery8_piecewise.py` | piecewise-profile fuzzer, with `expm` composed across segments as an **exact** oracle |
| `battery9_generic.py` | the generic user-Hamiltonian entry points |
| `battery10_coverage.py` | the whole Earth/solar surface a user actually touches, 164 configurations |
| `gamma_slack_sweep.py` | the population behind `adiabatic.GAMMA_TO_ERROR` |
| `timing.py` | the alternating timing harness, with two controls the change cannot touch |
| `run_notebooks.py` | re-executes notebooks 02/03 and hashes every embedded figure |

`crosscheck_acceptance.py` deliberately does **not** import `cross_check_strategies`: that
function does not exist on the pre-fix tree, and an acceptance test that only runs where the fix
already is tests nothing. It reimplements the engine forcing with the `battery3.py` spy pattern
so it runs anywhere:

```bash
git worktree add /tmp/mainwt 978663a
PYTHONPATH=/tmp/mainwt/src python crosscheck_acceptance.py    # 7/7 known silent misses seen
python crosscheck_acceptance.py                               # regression check
```

`bitident.py` and `attribute.py` take an output path and are meant to be run twice:

```bash
git worktree add /tmp/mainwt main
PYTHONPATH=/tmp/mainwt/src python attribute.py attr_main.npz
python attribute.py attr_branch.npz
```

Every trap listed in `../HANDOVER_ADVERSARIAL_VALIDATION.md` is encoded once in `harness.py`
rather than left to be re-hit — in particular `vcc_func_from_rho_func`'s 7th positional
argument (`density_is_of_number_of_electrons`, **not** `nubar`), which is passed by keyword
here so the mistake is unrepresentable, and the one-parameter `H` closure, which is always
factory-built.
