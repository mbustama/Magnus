# Handover: the pre-publish audit, second session

**Branch** `pre-publish-audit`, HEAD `8c3b19e`, **nothing pushed**.
**Version is `1.0.0`** — the release tag must be `v1.0.0` exactly; `publish.yml` fails on any other.

This continues `HANDOVER_PREPUBLISH_AUDIT.md`, which is still accurate for everything it
describes. **Read that one first for §4 (accepted trade-offs) and §5 (how the audit went
wrong four times).** This file covers the five commits since, and it exists mainly so the
next session does not redo work.

Working tree is **dirty on purpose**: `notebooks/make_notebooks.py` and
`docs/source/averaged_probability.rst` are edited and held back until the notebook rebuild
(§3.1), so the generator and its regenerated outputs land in one commit.

---

## 1. Do not redo these

Everything in this section was done, measured, and is committed. The measurements are in the
commit messages; `git log` is the detailed record.

| | outcome |
|---|---|
| **Read all 27 notebooks' prose** (§3.1 of the old brief) | Done. Three findings, all fixed. Do not re-read. |
| **SEV_TOL** (§3.4, left unresolved) | **Resolved** — `O47`. See §2.1. |
| **PREM recipe in the notebooks** (§3.2) | Done, and it was **five** notebooks, not the three the old brief lists: 02, 03, 04, **05, 07**. |
| **`13-bitident` harness bug** (§1) | Fixed — `O48`. Verified by running the exact `run_job` invocation shape. |
| **Clean-room install** (§1, "worth one re-run") | Re-run: **1149 passed, zero failures**. N43 confirmed to have closed the 12 failures. Do not re-run. |
| **The `angles` feature** | Complete across 95 functions — `O49`. See §2.2. |
| **Feature-interaction sweeps** | Three sweeps, preserved in `docs/dev/interaction_sweeps/`. See §5. |

**Verified clean, no action needed, do not re-check:**

* Every documented `Default: X` matches its signature across seven modules (one known false
  positive: `2*pi` vs `6.28318…`).
* All twelve `S12…D31_BF_NUFIT_6_0` constants agree with `load_nufit_params` to 1e-9.
* No function anywhere in `src/magnus` takes a mixing angle without `angles` (AST scan).
* No shared helper hardcodes a function name in a message (the one that did is fixed).
* No `versionchanged` anywhere — the pre-1.0.0 convention holds.
* No `TODO`/`FIXME`/`XXX` in `src/`.
* `validate_input_battery` only type-checks `osc_params`, so it is convention-agnostic.

---

## 2. What the five commits did

### 2.1 `O47` — SEV_TOL: the two calibrations never disagreed

The old brief's §3.4 recorded contradictory evidence and recommended leaving it. **The
contradiction was an artefact.** The two calibrations used *different spectrum families under
the same "scale 1e2" label*: `[-s, -s(1-d), s]` spans `2s` and gives `m ≈ 0.44 s²`, while
`[0, d, S]` spans `S` and gives `m ≈ 0.11 S²` — a factor of four in `m`. Compared at equal
`m` they agree.

Second finding: **no value of the gate could have made the old 2e-13 claim true.** At
`m = 1111` — the cell `test_sev_tol_sits_inside_its_calibrated_window` pins as
must-stay-on-kernel — 11 of 1200 random bases exceed 2e-13. The gate must sit above that
cell, so lowering cannot rescue the claim.

`SEV_TOL` stays at `1e4`; **nothing about the kernel changed**. The docstring now states the
measured guarantee: 5e-12 across the admitted range, 5e-13 in the `m ≤ 1.1e3` corner, each
about twice the worst measured. `calibrate_sev_tol.py`'s header is now `RESOLVED`.

### 2.2 `O49` — the `angles` convention

`'sin'` (default), `'sin2'` (sine **squared**), `'rad'`, `'deg'`, on 95 functions across
`hamiltonians`, `oscprob`, `oscprobstd`, `magnus prob` and `globaldefs.load_nufit_params`.
Under `'deg'` the CP phases are degrees too; under the other three they stay radians.
Matches the keyword NuOscProbExact uses, so the two codes can be driven from one parameter
set. **The default is a pass-through, so nothing that does not ask for it changes.**

Verified: 16 wrappers × 4 conventions end to end, worst 2.2e-13. `'sin2'` and `'rad'`
round-trip exactly; `'deg'` cannot (degrees↔radians goes through π), and its 2.8e-17 input
error is amplified by accumulated phase to 1.6e-12 on a solar ray — nine orders below the
1e-3 default tolerance. `load_nufit_params` reproduces the published NuFit table:
`sin2` 0.3088/0.47/0.02248, `deg` 33.76/43.28/8.62 and `dCP` 212.

**Four silent-wrong paths were found and closed while building it**, none visible to any
existing gate — see §4.

Rode along in the same commit (same files, could not be split without leaving an
intermediate commit failing its own tree guard): the CLI could only reach the NuFit **6.0**
parameter sets although the default had moved to 6.1; seven messages in
`validate_input_battery` hardcoded `osc_prob_matter_std_potential` instead of using the
`source_func_name` they were handed; and `tests/test_routine_listings.py`.

### 2.3 `O50` — the sterile matter term and the Earth's density describe different media

L37 made the Earth's `Y_e` a function of radius and derived `r = (1 − Y_e)/Y_e` from it layer
by layer for the **density**. The sterile states' entry in the matter projector, `r/2`,
cannot follow — it is one matrix for the whole chord — so it kept taking the caller's scalar,
whose default is `1.0`: isoscalar, `Y_e = 0.5`, **precisely the composition L37 replaced**.

Measured at `costhz = -0.95`, `s14 = 0.15`, `s24 = 0.10`, `D41 = 1 eV²`: the isoscalar
projector differs from one built with the core's own `r = 1.1478` by **2.1e-02** in
`P(νμ→νμ)` — twenty times the default tolerance, flat in tolerance, and silent. Three
flavours never saw it (empty sterile block), which is why `A2b` survived a max-effort review.

**Reported, not resolved**: no single `r` is right for a chord crossing iron and rock. The
Earth wrappers raise `globaldefs.SterileMatterCompositionWarning` when the two disagree by
more than 2%, naming the **path-averaged** ratio for the chord asked for. **No number changes
unless the caller acts.** Twelve wrapper docstrings had also asserted the ratio "must match
the value given to `vcc_func_from_rho_func`", which on the Earth path is unsatisfiable.

*Design note:* a warning keyed on the *range* over the four layers is useless — the ocean's
`r = 0.80` drags it below the isoscalar 1.0, so the default sits inside it and is never
questioned, while on a core-crossing chord the ocean is 3 km of 12 000. Path-averaging
weights each layer by how much of the trajectory is in it and gives one number to pass.

### 2.4 `O51` — three keywords that were reachable and broken

* **The Sun could not be told what it is made of.** `osc_prob_{4,5}nu_sun` and their
  `_nsi`/`_liv` variants never exposed `ratio_number_neutrons_to_protons` and delegated
  without forwarding it — two passing a hardcoded `1.0`. The Sun is hydrogen-rich:
  `Y_e = (1+X)/2` runs ~0.68→0.88 and `r` ~0.47→0.14, so isoscalar `1.0` is **outside the
  physical range entirely**, unlike the Earth where it at least sits among the layers. All
  six now take and forward it: on `osc_prob_5nu_sun_liv`, moving to the Sun's own `r ≈ 0.29`
  changes the averaged survival probability by **4.5e-03**. Exposed rather than defaulted,
  because the solar profile is a fit to the electron *number* density — `Y_e` is already
  inside it and there is nothing to derive `r` from.
* **`t_breakpoints` was unusable on every Earth wrapper.** They place slab edges on the PREM
  crossings themselves, so a caller's argument collided in `**kwargs` and raised
  `got multiple values for keyword argument 't_breakpoints'` two layers down — while the
  package's own unrecognised-keyword message advertises it as forwardable. Now **merged**
  with the PREM crossings; `t_slab_edges` remains the way to place every edge yourself.
* **`n_jobs` is not a pure performance knob.** Splitting slabs across workers changes
  arithmetic order, and the refinement ladder's stopping test compares successive levels.
  Serial vs two workers on a 3ν PREM chord: **1.2e-03 apart at the default `rtol = 1e-3`**,
  6.6e-08 at 1e-6, 5.6e-11 at 1e-9. Within tolerance, not bitwise. Now documented.

**A test stopped me and was half right — read this before touching it.**
`test_the_sun_takes_no_electron_fraction` encoded L37's decision that the Sun exposes no
composition parameters. That conflates two jobs of `ratio_number_neutrons_to_protons`: the
average nucleon mass for a mass→number density conversion, which genuinely never runs on a
solar profile, and the sterile projector entry, which has nothing to do with it. L37 was
right about the three density parameters and right at 2ν/3ν; it went one step too far at
4ν/5ν. The test now says exactly that, with the measurement. **It was refined, not reversed.**

---

## 3. What is left

### 3.1 The notebook rebuild — **needs an idle machine**, blocks two files

**Nine notebooks must rebuild together.** `make_notebooks.py --only` refuses to finish while
any *other* notebook on disk is stale, so they cannot be done piecemeal, and its staleness
check runs *after* execution — a partial run writes its notebooks, exits 1, and skips
`extract_gallery()`, leaving the gallery PNGs stale.

```bash
docs/dev/run_capped.sh -m 6G -- python3 notebooks/make_notebooks.py \
    --only 02,03,04,05,07,13,14,23,25
```

(`--only` takes a **comma list** of fragments, not repeated flags.)

Why each is stale:

* **02, 03, 04, 05, 07** — `num_density_e_func_prem` now uses the layered PREM composition.
* **13** — section 4 rebuilt around `average=True` vs the adiabatic closed form (§3.2).
* **14, 23** — prose only, but the generator changed so the files are stale.
* **25** — sections 4/5 pin `YE_UNIFORM = 0.5` on all three sides.

**Timing-sensitive**: notebook 25 prints ms and µs/probability via `best_of`/`timed_batch`.
Notebook 25 executed in **38.5 s** in the overnight log; all 27 took 36 min.

Before: check `uptime` is quiet, kill any stray background jobs.
During: **do not interrupt** — an interrupted run leaves notebooks written-but-unexecuted,
which looks exactly like stripped outputs.
After:

* notebook 25 prints a `control ratio` from its interleaved control — check it is near 1.0;
* re-read nb25 §4/§5 prose for numbers that moved with the composition change: the
  `max |Magnus - NuOscProbExact|` lines, "NuOscProbExact sits at 2.6e-04, some six hundred
  times above that floor", the ~400x cost claim, and the §5 summary table;
* nb13's window-convergence cell should print 0.592430 / 0.592685 / 0.596496 / 0.602295, and
  its adiabatic table worst 3.33e-16;
* then `python3 tests/test_file_tree.py --write` if any file was added, the docs build, and
  commit `make_notebooks.py` + `docs/source/averaged_probability.rst` **together** with the
  regenerated notebooks.

### 3.2 Why notebook 13 was rebuilt (context for reading the diff)

Notebooks 13, 14, 23 and `averaged_probability.rst` all asserted a solar result of
"instantaneous 1.4e-03 → averaged 2.6e-05, a **53x** reduction". The script the docs cite as
their source (`avg_check.py`) now gives **8.889e-04 → 6.051e-04, 1x**. It went stale at
`0bf3a40`, which is **on `main`** — pre-existing, not audit-caused.

The deeper problem was not the numbers. Notebook 13's "averaged probability" was the **mean
of a 121-point scan over six *vacuum* oscillation lengths** on a trajectory in matter. That
estimator has no converged value: the truth's window mean drifts 0.592430 (6 L_osc) →
0.602295 (48), *away* from a limit, because a wider window also averages over changing
density. Meanwhile `average=True` matches the textbook adiabatic MSW closed form
(both ends in matter) to **3.33e-16** across 1–20 MeV.

So section 4 now shows the trap, proves the estimator unconverged, then gives the exact
route. `avg_check2.py` was re-run for the other three rows of the docs table; all four rows
were replaced with measured values and the unsupported "twentyfold" rule dropped.

### 3.3 Open decisions — **for the user, not the next session to decide alone**

1. **Should the Earth's sterile projector `r` default to the path-averaged value instead of
   `1.0`?** Physics. It would make the default self-consistent, but changes 3+1/3+2 Earth
   results by up to 2.1e-02 without the caller asking, on an already-audited branch. If it
   changes: notebook 07 and notebook 25 §5/§7 need re-execution, and §7's frozen dataset is
   built on the uniform convention so it would need `electron_fraction=0.5` pinned to stay a
   like-for-like comparison.
2. **`anim_earth.gif`** still shows the old uniform composition. `RENDER = True` in notebook
   27, ~9 min for that scene, then `tools/make_demo_video.py --shrink`. **Safe on a busy
   machine** — rendering is deterministic, not a measurement. It changes tracked files under
   `img/`, which makes notebook 27's compression table (`anim_earth.gif | 20.0 MB | 2.61 MB |
   7.6x`) stale, so 27 joins the rebuild set.
3. **Release**: push, PR, merge, tag **`v1.0.0`**. `D15` enforces tag == `pyproject` version.
   The 3.13 CI job (`C9`) has never run and is the early warning for numba lacking a wheel.
   Note `pyproject`'s classifiers **deliberately** stop at 3.12 — a classifier is a promise
   that something tested that version. Do not "fix" that; it is documented in `pyproject`.

---

## 4. The failure mode this session added

The old brief's §5 says its four self-inflicted failures were all *scripted edits to
structured text*, all caught by the docs build and none by the suite. This session found a
fifth category that **the docs build cannot see either**:

> **A parameter accepted and silently dropped.**

Four functions took `angles` and ignored it (`pmns_mixing_matrix`,
`hamiltonian_{3,4,5}nu_vacuum`). The Sun did it with the composition. `t_breakpoints` did the
inverse, refusing a keyword it advertised. **Tests passed, docs built, ruff was clean for all
of them.** The only thing that found them was an AST check for "declares the parameter, never
reads it", now in `tests/test_angles.py`.

Two rules that would have prevented the whole class:

1. **When adding a keyword, check three texts agree** — signature, docstring, *and body*. The
   third is the one no existing gate covers.
2. **Cross the axes.** Every defect of consequence this session found was an interaction:
   layered `Y_e` × sterile flavour count, cache key × convention, defaults helper ×
   convention, Earth wrapper × `t_breakpoints`.

Also learned, twice: **when a sweep disagrees, suspect the harness first.** `strategy` in
vacuum and bitwise `n_jobs` agreement were both my expectations being wrong, not the library.

---

## 5. New guards and tools

**In the suite (run in CI):**

* `tests/test_angles.py` — 34 tests: the four conventions agree at builders, wrappers and
  closed forms; every guard fires and stays silent when it should; `load_nufit_params`
  reproduces the published table; the mixed-convention defaults bug is pinned *with its
  size* (0.143 in probability); the sterile-composition warning; the solar ratio; the
  `t_breakpoints` merge. Includes the AST check for the "accepted and dropped" class.
* `tests/test_routine_listings.py` — 21 tests: every public function is named in its module's
  `Routine listings`. Found three omissions in `magnus.magnus` the first time it ran.
  `oscprob` and `plotting` are pinned as deliberate exemptions so a section cannot silently
  disappear either.

**In `docs/dev/interaction_sweeps/`** (kept, not in CI — see its `README.md`): the three
cross-feature sweeps and two structural checkers. Reuse the grid rather than re-deriving it.

---

## 6. Quick reference

```bash
# fast gates (~1 min each)
python3 -m ruff check src/ tests/ docs/ notebooks/make_notebooks.py
python3 -m pytest tests/test_file_tree.py -q
python3 docs/regen_cli_help.py --check        # regenerate after ANY cli.py change

# the gate that catches scripted-edit damage (~2 min)
cd docs && rm -rf build source/api && python3 -m sphinx -b html -n -W --keep-going source build

# full suite, ~10 min -- currently 1205 passed
python3 -m pytest tests/ -q -n auto

# regenerate the file trees in README.md and installation.rst
python3 tests/test_file_tree.py --write        # NOT `pytest --write`, which is rejected
```

**Traps that cost real time, this session and last:**

* `make_notebooks.py` **executes** notebooks; interrupting leaves them written-but-unexecuted.
* `git add -A` while a notebook build runs stages half-written notebooks.
* `overnight_audit.sh` **refuses a dirty tree** — correct behaviour, but it means the script
  is blocked until you commit.
* A new test file must be **git-tracked** before `test_file_tree.py` passes (it asserts
  listed == tracked). `git add` it before running the suite.
* Adding a public name to `globaldefs` does **not** document it — that module has an explicit
  `__all__` allowlist which autoapi respects. Two names rendered nowhere until they were
  added to it, and the docs build stayed green throughout. **Grep the built HTML** for any
  new public symbol.
* A private module (`_angles.py`) is never rendered, so cross-references inside it are never
  checked by `-n -W`.
* `angles` must go **last** in a signature, and **before** `**kwargs` where one exists —
  appending blindly produced 60 syntax errors in one pass.
