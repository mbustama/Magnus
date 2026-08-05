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

Added by the physical-profile programme (`../HANDOVER_PHYSICAL_PROFILES.md`), which asks whether
any of the earlier findings reaches a real user:

| file | what it does |
|---|---|
| `physical_profiles.py` | the physically-motivated population: tabulated profiles with interpolation kinks, the real BS05(AGS,OP) solar model, a supernova shock front, Kolmogorov turbulence, and Earth with a non-PREM crust. Every family carries its own trajectory, energy band and a `provenance` string saying how physical it actually is |
| `validate_physical.py` | checks each profile has the shape it claims — jump factors against the literature formulas, spectral index, kink placement, resonance on the trajectory — *before* anything measures with it |
| `physical_battery.py` | the two questions no existing script covers: `sub_grid` (P2/P3, `find_hidden_features` rates) and `seam_cost` (P4, cumulative-vs-hybrid cost on physical profiles) |
| `bs05_agsop.dat` | the BS2005-AGS,OP standard solar model table, from `sns.ias.edu/~jnb/SNdata/`. **Third-party data** -- see `SOURCES.md` for provenance, citation and terms |
| `attribute_physical.py` | attribution for the silent misses the physical population produced: which engine answered, whether it is a knife-edge, and whether `t_breakpoints` or a tighter request cures it |
| `bs05_energy_band.py` | whether the BS05 silent miss reaches a real solar-neutrino energy, with the oracle verified at each one |
| `avg_check.py` / `avg_check2.py` | whether an error survives PHASE AVERAGING -- the observable for solar and supernova physics. Collapses the solar error 53x and the turbulence error 23x, and leaves the shock error untouched |
| `shock_silent_band.py` | maps shock width x energy for "wrong, unflagged AND unwarned". 18 configurations, 2 outside tolerance, both warned, none silent |
| `alias_cost.py` | what an aliasing check would COST, measured before it was written |
| `alias_fp.py` | what it would FIRE on -- 44 of 45 realistic scans, which is why it ships as a `strategy_info` statistic and not a warning |

`fallback_quality.py`, `warn_fp.py`, `resolution_fp.py` and `crosscheck_acceptance.py` all take
a `--physical` flag that swaps the population and leaves the measurement logic alone. They write
to a separate `*_physical.npy` so a physical run and a synthetic one cannot overwrite each
other's rows:

```bash
python validate_physical.py            # first: are the profiles what they claim?
python warn_fp.py --physical           # P1, the headline: any silent miss?
python physical_battery.py sub_grid    # P2/P3
python physical_battery.py seam_cost   # P4
```

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
