# Handover: expanding the Magνs notebook set

**Written:** 2026-08-08, at the close of a long session. **Read §§0–2 before touching anything.**

**The task for the next session.** Look at NuOscProbExact's twenty example notebooks, report
which could be adapted for Magνs — *adapted*, not ported, since the two packages solve different
problems — and propose additional notebooks that only Magνs can support. **Produce the report
first and wait for the user's decision before writing any notebook.** The user asked for
suggestions of your own as well as adaptations.

---

## 0. Verify the base before starting

```bash
git -C ~/Research/magnus branch --show-current      # notebooks
git -C ~/Research/magnus log --oneline -1           # "Hand over the notebook-expansion task"
git -C ~/Research/magnus status --porcelain         # empty
python -c "import sys; sys.path.insert(0,'src'); import magnus.oscprob as o, magnus.magnus as m; \
  print(o.MIN_EFFECTIVE_REFINEMENT, o.BATCH_WORKING_ENTRIES, m.USE_PALINDROME)"
# must print: 1.25 65536 True
```

**Baseline:** 894 tests pass (~20 min). `ruff check src/magnus/ tests/ docs/
notebooks/make_notebooks.py notebooks/make_shock_reference.py` clean. `cd docs && make html
SPHINXOPTS="-n -W --keep-going"` clean **after `make clean`**. `python docs/regen_cli_help.py
--check` in sync. `python tests/test_file_tree.py` says up to date.

**Branch `notebooks` is ahead of `main` and fully pushed** as of this brief. `main` does not
carry any of it yet; the branch has not been opened as a pull request.

---

## 1. THE CRITICAL FACT: the notebooks are generated

**`notebooks/make_notebooks.py` (6600+ lines) builds all fourteen `.ipynb`. Edit the generator,
never the notebook.** Anything written into a `.ipynb` by hand is lost on the next build.

```bash
python notebooks/make_notebooks.py                 # build + execute + store outputs (~37 min)
python notebooks/make_notebooks.py --no-execute    # build only -- WIPES ALL STORED OUTPUTS
```

**`--no-execute` is destructive.** It rewrites every notebook from source with no outputs. If you
use it and then execute only the one you changed, the other thirteen are left blank and
`test_file_tree`-adjacent CI (`.github/workflows/notebooks.yml`) will fail on the
stored-outputs check. Recovery: `git checkout <commit> -- notebooks/NN_*.ipynb` for the ones you
did not change, then verify their cell sources still match the generator. That happened this
session and cost a re-run.

Structure of the generator: `md()`, `code()`, `notebook(title, intro, cells)`, a `books` dict,
`READING_ORDER` (filename, title, blurb), `add_footers()`, `extract_gallery()`, `build()`.

* `add_footers()` **asserts `READING_ORDER` and `books` have the same keys**, so a new notebook
  fails the build until it is added to both.
* `extract_gallery()` lifts figures out of the executed notebooks into `img/gallery/`, keyed by
  `(notebook, index of the PNG output within it)`. Adding a notebook before an existing gallery
  entry does not shift anything, but reordering *cells* within a gallery source notebook does.
* Adding a notebook also needs: an entry in `docs/source/tutorials.rst` (grouped by purpose, not
  numbered), and a line in `tests/test_file_tree.py`'s `TREE` — the suite fails otherwise.

### Execution cost, measured 2026-08-08

Total **~37 min**. Per notebook: 14 → 330 s (was 933 before its oracle was frozen), 12 → 411 s,
03 → 318 s, 02 → 289 s, 07 → 218 s, 13 → 197 s, 06 → 150 s, 11 → 147 s, 04 → 64 s, 01 → 60 s,
08 → 25 s, 10 → 11 s, 05 → 12 s, 09 → 7 s.

**Any new notebook adds to this and to every CI run.** `notebooks.yml` is paths-filtered
(`notebooks/**`, `src/magnus/**`, `pyproject.toml`) so docs-only changes skip it, but a notebook
change runs the whole set. Budget accordingly: a notebook that takes ten minutes needs to earn it.

### The frozen-oracle pattern, for anything expensive

Notebook 14 spent ~600 s recomputing a `solve_ivp` ground truth that **cannot change** — it
depends on the profile and the energy, not on Magνs. It is now computed once by
`notebooks/make_shock_reference.py` into `notebooks/shock_reference.json` as hexadecimal floats
(exact round-trip), and the notebook loads it. **All twelve numeric outputs were verified
identical to the live-oracle run before the change was accepted.**

Two rules that came out of building it, both paid for:

1. **Do not transcribe the physics into the reference generator.** The first version hand-copied
   the shock profile and got it wrong — a power-law rarefaction instead of the Fogli et al. form,
   and the wrong density normalisation. `solve_ivp` overflowed, which is why it was caught; had
   the numbers merely been *plausible* it would have frozen an oracle for a profile the notebook
   does not use, and every comparison against it would still have looked fine. The generator now
   **executes the notebook's own cells** out of `make_notebooks.py` up to the point the reference
   is needed. One definition of the physics, read rather than copied.
2. **Store a fingerprint.** The file carries the electron density sampled along the ray; the
   loader raises if the profile it just built does not match. A stale oracle that silently
   outlives a change to the physics is worse than no oracle.

**Do not freeze an oracle whose *timing* is part of the point.** Notebook 12 is the second-slowest
at 411 s and was deliberately left alone: its `run_case` times `solve_ivp` as part of the
comparison it exists to make, so freezing would replace a live measurement with a stale claim.

---

## 2. NuOscProbExact's twenty notebooks, and how they map

Source: `/home/mbustamante/NuOscProb/NuOscProbExact/notebooks/`, generated by that project's own
`make_notebooks.py` (5979 lines). Titles and blurbs from its `READING_ORDER`.

| # | Their notebook | Blurb | Magνs already has? |
|---|---|---|---|
| 01 | Basics | units, one probability, why to pass arrays | **yes** — 01 |
| 02 | Oscillations in vacuum | against baseline and against energy | **yes** — 02/03 |
| 03 | Matter, NSI, and LIV | constant density and two kinds of new physics | **yes** — 02/03/08/09 |
| 04 | Oscillograms | a 2-D map in one call | **yes** — 06 |
| 05 | Bi-probability plots | CP violation as an ellipse | **yes** — 05 |
| 06 | The Earth: PREM, chords, slabs | how a varying profile becomes exact pieces | partial |
| 07 | Probabilities through the Earth | zenith scans, Earth oscillogram, real baselines | **yes** — 04/06 |
| 08 | Unusual density profiles | castle walls, why the arrangement matters | partial — inside 02/03 |
| 09 | Performance | looping versus broadcasting, measured live | **no** |
| 10 | The paper's figures | reproduces arXiv:1904.12391 | n/a — no Magνs paper yet |
| 11 | Exact versus the approximations | where the familiar formulas break | **no** |
| 12 | Mass ordering and the octant | two open questions, how they show up | **no** |
| 13 | Antineutrinos, done properly | conjugate *and* flip; two ways to get it half right | **no** |
| 14 | Solar neutrinos and the MSW resonance | the adiabatic resonance, limits of slabbing | partial — 12/13 |
| 15 | Numerical edge cases | degeneracies; what returns a number instead of NaN | **no** |
| 16 | Four neutrinos, and a sterile state | 3+1 through SU(4), why the method stops | **yes** — 07 |
| 17 | Cross-checks with other codes | nuSQuIDS, Zaglauer–Schwarzer | **no** |
| 18 | The evolution operator and SU(n) coefficients | the machinery underneath | partial — 11 |
| 19 | Animated scenes | four sweeps as stills | **no** |
| 20 | An arbitrary Hamiltonian, through three profiles | a long-range force through three bodies | **no** |

Magνs's fourteen, for comparison: 01 introduction, 02 2ν, 03 3ν, 04 long baseline,
05 biprobability, 06 oscillograms, 07 sterile, 08 NSI, 09 LIV, 10 averaged, 11 matrix
exponential, 12 strategy, 13 tabulated solar, 14 supernova shock.

### My reading of what is worth adapting

**Strongest candidates**, in order:

1. **13 Antineutrinos, done properly.** Magνs has *no* antineutrino notebook, and this is the
   single most defect-prone convention in the package — the sign flip was applied **twice** at one
   point, giving antineutrinos a positive matter potential and plausible-looking wrong answers.
   The convention is now written down (`docs/source/methodology.rst`, `Conventions` section, anchor
   `_conventions`) but never *demonstrated*. Their framing — "conjugate *and* flip, and two ways to
   get it half right" — transfers directly and is worth more here than there.
2. **11 Exact versus the approximations.** Magνs's analogue is stronger: it can show where the
   *two-flavour closed form* and the constant-density formula break against a real varying profile,
   which is exactly what the package is for. Reuse `oscprobstd`'s closed forms as the "familiar
   formula" side.
3. **09 Performance.** Magνs has a lot to show that they do not: the palindrome (1.4–1.67× on an
   expensive `H_func`), array-vs-loop (4.6× measured), the batched separable engine, and
   `BATCH_WORKING_ENTRIES` (1.19–1.38×). All measured this session; see `docs/source/index.rst`'s
   Performance section for the numbers and their populations.
4. **15 Numerical edge cases.** Degenerate spectra, `dCP` at 0 and π, zero baseline, a single
   slab, `rtol=atol=None`. Magνs-specific additions: what `ToleranceNotAchievedWarning` means, what
   `convergence_info['tolerance_achieved']` reports, and the case where two levels agree
   coincidentally (the **0.855** error on a sawtooth density that `strict_convergence` exists for).
5. **12 Mass ordering and the octant.** Straightforward and physically useful; Magνs has both
   NuFIT parameter sets (`OSC_PARAMS_NU_FIT_6_0_SK_NO` / `..._SK_IO`).
6. **08 Unusual density profiles.** Magνs already touches castle-wall and noisy profiles inside 02
   and 03, but a dedicated notebook could carry the `t_breakpoints` story properly — including the
   measured asymmetry that on a **scan** breakpoints are an established cure while on a **single
   point** they improved 7 of 18 configurations and *worsened 11*.

**Weak or not applicable:** 10 (no Magνs paper yet — ask the user), 16 (Magνs's 07 already covers
4ν and 5ν, and Magνs has no SU(4) ceiling to explain), 17 (would need an external code installed;
Magνs's cross-checks are against `solve_ivp`, already in 12/13/14), 19 (animations are expensive
and the CI budget is already 37 min), 06 and 18 (largely covered by Magνs's 11 and by
`docs/source/methodology.rst`).

### Notebooks only Magνs can support — my own suggestions

These have no NuOscProbExact counterpart because they are about *this* package's machinery:

* **"Which engine answered, and why."** `strategy_info` reports the engine; Magνs has five
  (hybrid, ip_exp, separable, cumulative, general ladder) with a documented dispatch order and
  measured reasons for it. `cross_check_strategies` runs several and reports their spread —
  and its own docstring example was **broken** until this session (it passed `s12` to a
  two-flavour call, so nothing ran, hidden behind `# doctest: +SKIP`).
* **"What `rtol`/`atol` actually promise."** The most valuable single notebook I can think of.
  It is a *stopping criterion*, not an error bound; usually conservative; occasionally not. Show
  `convergence_info`, show `n_slab_edges` versus `n_slabs` (a nominal 2→3 slab step is a 16→17
  **edge** step on an Earth chord), and show the defect fixed in PR #35 where the ladder certified
  an agreement between two nearly identical grids. See `docs/source/implementation_details.rst`,
  anchor `_what-rtol-atol-control`.
* **"When averaging rescues you and when it does not."** Magνs's 13 and 14 are already the two
  halves of this; a short notebook could put them side by side. Error that is a *phase* falls 53×
  under averaging; error that is an *envelope* does not move.
* **"Bring your own Hamiltonian."** `osc_prob_earth`/`osc_prob_sun` take a user `H_func`. Cover
  the vectorisation trick (`[..., None, None]`, worth 4.6×), `ScalarHamiltonianWarning`, and the
  palindrome declaration that the Earth entry points make on the caller's behalf.
* **"Reproducing a published figure."** Ask the user whether there is a paper or a figure they
  want reproduced; this is their call, not mine.

---

## 3. Traps from this session, all paid for once already

* **Run every docstring/notebook example before writing it down.** Two of fifteen candidate
  examples were wrong when written from the signature (`oscprobstd.osc_prob_2nu_*_std` return a
  2×2 matrix, not a tuple). Worse, a *shipped* example ran nothing at all behind `+SKIP`.
* **`--no-execute` wipes stored outputs of all fourteen.** See §1.
* **Do not hand-transcribe physics between files.** See §1.
* **A `:doc:`/`:func:` reference that resolves to nothing renders as plain text and Sphinx does
  not report it** unless run with `-n`. CI now builds with `-n -W`, and `conf.py` has
  `nitpick_ignore_regex` for the ~2470 numpydoc type-string false positives that would otherwise
  make `-n` unusable. **If you add a notebook that references a new symbol, `-n` will catch a
  typo; do not disable it.**
* **A public symbol missing from `__all__` is invisible to sphinx-autoapi**, which silently
  breaks every reference to it. Twenty-one were found this way: eighteen tuned constants,
  `PhaseAveragingWarning`, and the `Term`/`Word` type aliases.
* **A cell that exactly fills its column breaks an RST simple table** ("text in column margin").
  Establish it by running `docutils.core.publish_doctree` on the table, not by eye.
* **The file tree is generated.** `tests/test_file_tree.py` — add new files to `TREE` or the suite
  fails. `python tests/test_file_tree.py --write` regenerates both README.md and
  installation.rst. Its writer is destructive if block detection is wrong: the first version
  matched ``` exactly, missed that the fence is ```` ```text ````, and rewrote the wrong region.
* **`git stash` while a background pytest run is in flight invalidates it.** Re-run.
* **`gh` fails on an unreadable `/etc/gitconfig`**; `GIT_CONFIG_NOSYSTEM=1 gh ...` works.
* **Kill background waiters by PID.** `pgrep -f` matches the shell issuing it. Never leave two
  waiters on one condition — and never wait on a sentinel that a killed job will never print.
* **The user asks for brevity.** Lead with the answer; keep the measurements exhaustive and the
  prose short. See `~/.claude/.../memory/magnus-be-brief.md`.

---

## 4. Facts a notebook author will want, with sources

* **Conventions** — `docs/source/methodology.rst`, anchor `_conventions`. `P[nu_i][nu_f]`,
  initial flavour first. Antineutrino sign applied once inside
  `matter.vcc_func_from_rho_func`. Mass ordering carried by the **sign** of `D31`. Angles are
  `sin θ`, not `sin² θ`. Two flavours take `sth`/`Dm2`, not `s12`/`D21` — and passing the wrong
  names is **silently ignored**, not an error.
* **Units** — `docs/source/quickstart.rst`, anchor `_units-table`.
* **Recipes, all executed at docs build** — `docs/source/recipes.rst`. Eleven short ones; the
  fastest way to see a working call for anything.
* **Engines, dispatch and every tuned constant with its measured population** —
  `docs/source/implementation_details.rst`.
* **The palindrome** — worth 1.4–1.67× on an expensive `H_func` across an Earth chord, ~0.9× on
  plain PREM, nothing on a standard PREM *scan* (the separable engine already shares the profile).
  `magnus.magnus.USE_PALINDROME` disarms it. See `docs/dev/PLAN_PALINDROMIC_PROFILES.md` §3f.
* **Warning classes** — `ToleranceNotAchievedWarning`, `MagnusConvergenceWarning`,
  `UnmarkedDiscontinuityWarning`, `HybridCertificationWarning`, `PhaseAveragingWarning`. Measured
  false-alarm rates are in `implementation_details.rst`; `MagnusConvergenceWarning` is 76% false
  and is a statement about slab width, not about the answer.

---

## 5. Still outstanding, unrelated to notebooks

* **The physical-profile population has never been run against the palindrome gate** from PR #33.
  §6 of `PLAN_PALINDROMIC_PROFILES.md` asks for confirmation that the mirror declines on solar,
  supernova and shock profiles. There is a unit test for the solar case, not the population. **This
  is the only open item with correctness rather than documentation consequences.**
* **P4** — `docs/dev/adversarial_batteries/RUN_P4.md`, unrun since the robustness programme.
* **A spawned task** on `simpson`'s even-`n_tpts` parity: scipy applies an asymmetric
  last-interval correction when the interval count is odd, and Magνs's default
  `n_tpts_per_slab` is 50 — that parity.
* **arXiv / Zenodo badges** await real identifiers from the user.
* **Example coverage is 39%** (76 of 192 public entities). The remaining `hamiltonians*`
  variants (`_td`, energy-dependent) and ~54 `oscprob` wrappers were deliberately left: they
  restate what is already shown.
