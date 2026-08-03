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
