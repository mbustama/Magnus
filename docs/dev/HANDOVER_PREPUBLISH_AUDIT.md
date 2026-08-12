# Handover: the pre-publish audit

**Branch** `pre-publish-audit`, 47 commits, HEAD `541af1b`, working tree clean, **nothing pushed**.
**Version is now `1.0.0`** (not rc1), so the release tag must be `v1.0.0` — `publish.yml` fails on any other.

Read this before running anything. Several jobs here take hours and have already been run.

---

## 1. Do not re-run these

An overnight verification completed at **2026-08-12 07:40**, logs in
`~/magnus_overnight/20260812-040953/` (`SUMMARY.txt` first). All of it ran at `5ff3ee0`, three
commits behind HEAD; the three since are `N43`, `N44` and a docs fix, none of which touch physics.

| job | result | re-run? |
|---|---|---|
| **P4 seam-cost** (`physical_battery.py seam_cost`) | **OK, 153 min.** Cumulative cheaper than hybrid everywhere; worst ratio 0.21x at N=8, median 0.02x; "the 25 -> 8 change holds" | **No.** First completion ever; supersedes the discarded contaminated attempt |
| **all 27 notebooks** (`make_notebooks.py`) | OK, 36 min | Only if you change `make_notebooks.py` or `src/` |
| clean-room install | 12 failures — **cause fixed in N43** | Worth one re-run to confirm N43 closed it (~20 min) |
| coverage | 93.26% | No |
| PREM batteries 11, 12 | OK | No |
| docs `-n -W` | 1 warning — **fixed in N44** | Cheap (~2 min); it is the gate that catches most of my mistakes |

`overnight_audit.sh` re-runs the lot. **Do not** invoke it casually — it is ~3.5 h and phase 2
needs an idle machine. Individual jobs can be lifted out of the script.

### Known harness bug, not a code bug
`13-bitident` failed with `IndexError: sys.argv[1]`. `bitident.py` needs an output filename and
`overnight_audit.sh` calls it bare. One-line fix in the script.

---

## 2. What was done

47 commits, grouped by the audit item they close. Every commit message carries the measurement
behind it; `git log` is the detailed record and is worth reading before touching any of this.

**Two wrong-answer physics bugs, both in published entry points:**

* `A2b` — the NSI route built its standard matter term as a literal `diag([1,0,0,0])`, so sterile
  states carried **zero** where they carry `-V_NC = (r/2)V_CC`. With every coupling zeroed it must
  reproduce the standard route and did not: **5.2e-02 at 4nu, 5.1e-02 at 5nu**. Three flavours
  agreed all along, which is why nothing caught it.
* `A2c` — `osc_prob_Nnu_sun_liv` read the Sun's **electron number** density as a **mass** density,
  because it forwarded a flag defaulting to False. Off by **up to 0.69** in probability.

Both found the same way: **ask what happens when the new physics is switched off.** Neither is
reachable by testing NSI or LIV. `A2d` now sweeps all nine BSM families x four flavour counts (36
cases) so the *shape* is guarded, not just the two instances.

**Two deliberate, results-changing defaults** (both in the changelog with reproduction instructions):

* `K35` — the implicit oscillation parameters were NuFIT 6.0 while `load_nufit_params()` returned
  6.1. Unified on 6.1, **derived from the loader** so they cannot drift again. 4.0e-03 in probability.
* `L37` — the Earth's `Y_e` is now resolved per PREM layer (core 0.4656 iron / mantle 0.4957
  peridotite / crust 0.4952 granitic / ocean 0.5551 seawater), and the neutron-to-proton ratio is
  *derived* from it, `r = (1-Y_e)/Y_e`. **Up to a factor of four on core-crossing chords.**
  `electron_fraction=0.5` reproduces the old uniform composition.

**The Sun deliberately has no `Y_e` parameter** and this is not an oversight — its profile is the
standard exponential fit to the *electron number* density, so the mass-density conversion never
runs and `Y_e` is already inside the fit. `L37` removed four parameters from `sun_liv` that were
accepted and silently ignored. Physically `Y_e` does vary through the Sun, more than through the
Earth — `(1+X)/2`, so ~0.67 in the depleted core to ~0.87 in the envelope — and continuously, so
the Earth's layered treatment would be the wrong model even if the conversion did run.

**Also:** validation guards (`A1` non-positive energy silently returned the *antineutrino*
probability; `B3`/`B4` slab counts; `B5` NaN density mis-reported as a units error); release
metadata (`C6`–`C12`: classifier, `py.typed`, CI matrix, badges, `CITATION.cff`, licence headers);
publish gating (`D13`–`D15`: tests, `twine check`, **tag-vs-version match**); changelog (`E16`–`E18`,
`J34`, `K36`, `L39`); docs accuracy (`F19`–`F24`, `J32`, `J33`); tests (`G25`, `G26`, `A2d`);
`H28`/`I30` error-message convention; `I31` exported `matter_potential_projector`; `N43` made numba
required and `expm_backend` discoverable.

**Four commits repair damage this audit itself caused** — `I29`, `I30`, `L38`, `N44`. See §5.

---

## 3. What is missing

### 3.1 Read the notebooks (not just execute them) — **highest value**
`04-notebooks-all` re-executed all 27, so their *numbers* are current. **Nothing has read their
prose.** Two defaults moved results; a notebook can execute cleanly while its text asserts
something its new figures contradict. The six that call Earth wrappers are 09, 12, 19, 24, 25, 27.
This needs eyes, not a script.

### 3.2 Notebooks 02, 03, 04 teach the old composition
They define `num_density_e_func_prem` with `electron_fraction=0.5` hardcoded. They do **not** call
the Earth wrappers, so nothing disagrees today — but a reader following that recipe and then
comparing against `osc_prob_3nu_earth` gets an unexplained mismatch. Same shape as the notebook-12
reference `A2e` fixed, latent rather than active.

### 3.3 `anim_earth.gif` shows the old composition
Committed in `img/`. Re-render is `RENDER = True` in notebook 27, ~9 min for that scene, then
`tools/make_demo_video.py --shrink`. Raw renders go to `img/raw/` (gitignored); the shrunk copy in
`img/` is what is tracked.

### 3.4 `SEV_TOL` — **unresolved, and I left contradictory evidence**
Read the `UNRESOLVED` section in `docs/dev/calibrate_sev_tol.py` first. In short:

* my grid finds cells the 1e4 gate **admits** reaching 5.1e-13 absolute against `eigh`'s 6e-15 —
  past the 2e-13 budget the `SEV_TOL` docstring claims — smallest offender at m = 4.4e3;
* `HANDOVER_OVERHEAD.md` records the original calibration finding the first unsafe cell at 1.1e5;
* `test_expm_backend.py::test_sev_tol_sits_inside_its_calibrated_window` encodes the original and
  **fails if the gate is lowered to 1e3** — so 1e3 would decline cells the original measured safe.

I recommended lowering it, then found that test and reverted. Lowering is not a conservative tweak
if it contradicts a standing calibration; it is picking a side. **Mitigating fact:** instrumented
across a PREM chord, a constant-density call, a 60-energy Earth scan and a solar profile, the
severity actually reached is **m <~ 10** — a Magnus slab has `||Omega|| <~ pi` by construction, so
neither candidate gate ever fires in ordinary use. This is about whether a documented guarantee is
true, not about numbers users are getting. Deciding to leave it and correct the docstring instead is
a defensible outcome.

### 3.5 One intermittent test
`test_magnus_expansion.py::test_cached_eval_mode_bounds_its_per_interval_dictionary` failed once
under `-n auto --cov` and **passes standalone and under `--cov`** — verified both. Order- or
state-dependent under xdist, not a defect.

### 3.6 Release
Push, open the PR, merge, then tag **`v1.0.0`**. `D15` enforces tag == `pyproject` version, so any
other tag fails the publish job by design. The 3.13 CI job (`C9`) has never run — it is also the
early warning for numba lacking a wheel on a new interpreter, now that `N43` made numba required.

---

## 4. Accepted trade-offs (do not "fix" these)

* **`git clone` is 225 MB** (not the 709 MB the local `.git` reports — that is unpacked loose
  objects; `git gc` shrinks it locally and changes nothing for anyone else). History rewrite was
  considered and **declined**. `pip install` is 0.32 MB, ~700x lighter; nothing under `img/` reaches
  the wheel or sdist, verified by building and inspecting both.
* **numba is required** (`N43`). Runtime opt-out is unaffected: `magnus.magnus.EXPM_BACKEND='eigh'`
  or `expm_backend='eigh'` per call, both agreeing with the compiled path to 8.88e-16. The cost is
  that numba lags new interpreters, so a Python it has no wheel for makes the package
  uninstallable rather than merely slower.
* **PREM's ocean default** (`Y_E_OCEAN_PREM`) encodes a global-average ocean a land baseline does
  not cross. Documented, not guessed; pass `electron_fraction_ocean=Y_E_CRUST_PREM` for one.
* **Crust vs mantle differ by 0.1%** — the crust parameter exists for explicitness, not effect.

---

## 5. How this audit went wrong, four times

Worth reading before making a similar change, because the mechanism repeated:

`I29` — F24's docstring insertion sliced the Parameters section at the first *word* "Returns",
which occurs inside parameter prose, corrupting 24 docstrings.
`I30` — H28 used `ast.walk` (breadth-first) to find "the first string" in a concatenation and put
the prefix mid-message in 12 places. f-strings needed separate handling again.
`L38` — L37's docstrings went in through a plain Python string, so `\rm`, `\frac` and `\rangle`
became control characters. Python emitted `SyntaxWarning: invalid escape sequence` and I read past it.
`N44` — N43 added a `:func:` reference to a private name, which autoapi does not document.

Every one was a **scripted edit to structured text**, and every one was caught by the **docs build**,
never by the suite — none of it is executable. Two rules that would have prevented all four:

1. **Raw strings** when writing LaTeX or RST through a script, and check for stray control bytes.
2. **Run `make html SPHINXOPTS="-n -W --keep-going"` before committing**, not after. `lint.yml`
   runs it on every push, so CI would have caught these — but only after a push.

---

## 6. Quick reference

```bash
# fast gates (~1 min each)
python3 -m ruff check src/magnus/ tests/ docs/ notebooks/make_notebooks.py
python3 -m pytest tests/test_file_tree.py -q
python3 docs/regen_cli_help.py --check

# the gate that catches scripted-edit damage (~2 min)
cd docs && rm -rf build source/api && python3 -m sphinx -b html -n -W --keep-going source build

# full suite, ~9 min
docs/dev/run_capped.sh -m 6G -- python3 -m pytest tests/ -q -n auto

# rebuild only the notebooks you changed -- NEVER a bare make_notebooks.py unless you mean all 27
docs/dev/run_capped.sh -m 6G -- python3 notebooks/make_notebooks.py --only 12_magnus_adiabatic_hybrid_strategy.ipynb
```

**Traps that cost real time here:**
`make_notebooks.py` **executes** notebooks; interrupting it leaves them written-but-unexecuted,
which looks exactly like "regeneration stripped the outputs". `--only` also verifies every *other*
notebook still matches the generator. Adding a notebook changes its neighbour's footer, because
footers are built from `READING_ORDER`. `tests/test_file_tree.py --write` regenerates the README and
installation trees; pytest alone only compares. And `git add -A` while a notebook build is running
stages half-written notebooks — nearly committed that.
