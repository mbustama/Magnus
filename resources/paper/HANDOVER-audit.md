# Handover: run the audit on the Magνs CPC paper

**Your job is one thing: run the audit in
`resources/paper/audit-criteria.md` over `resources/paper/main.tex`, to a clean
report.** Everything else described here is already done. Do not rebuild it.

Written 2026-08-30, at the end of the session that wrote the paper. Figure 2 was
rebuilt after the rest of this brief was drafted; §7 records what changed and is
the one place where the paper moved after the audit criteria were last swept.

---

## 1. What exists, and where

| What | Path |
|---|---|
| The paper | `~/Research/magnus/resources/paper/main.tex` |
| The audit you are running | `~/Research/magnus/resources/paper/audit-criteria.md` |
| How to build it, and every convention | `~/Research/magnus/resources/paper/README.md` |
| Bibliography | `~/Research/magnus/resources/paper/refs.bib` (199 entries) |
| Figures | `~/Research/magnus/resources/paper/figs/` (nine PDFs, tracked past a global `*.pdf` ignore) |
| The notebook that produces every figure | `~/Research/magnus/notebooks/28_magnus_paper_figures.ipynb` |
| **Which is generated from** | `~/Research/magnus/notebooks/make_notebooks.py` — **edit the generator, never the `.ipynb`** |
| The sibling paper, for style and shared references | `~/Research/NuOscProb/NuOscProbExact/resources/paper/main.tex` |

Build: `pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex`.
Four passes. The paper currently compiles at **23 pages, zero undefined
references, one overfull hbox**.

Regenerate all figures: `python notebooks/make_notebooks.py --only 28_magnus_paper`.
Roughly three minutes. It executes the notebook and stores its outputs.

---

## 2. The author's preferences, measured rather than assumed

All of this is already applied to the draft. Keep it when you edit.

**Writing style**, measured over 57 000 words of his arXiv papers *and* over his own
CPC code paper. The genre matters: his physics-results letters run at 8.6 first-person
per 1000 words, but `NuOscProbExact` — same author, same genre as this — runs at 3.6.
**Use 3.6, not 8.6.** Against that target the current draft sits at:

| marker, per 1000 words | this paper | `NuOscProbExact` |
|---|---|---|
| first person | 2.9 | 3.6 |
| em dashes | 6.5 (was 9.9) | 4.6 |
| parentheses | 9.4 | 16.6 |
| `\ie` / `\eg` | 0.1 | 0.6 |

He reaches for parentheses where I reached for em dashes. Converting an aside from
dashes to parentheses fixes both numbers at once. Sentences median 18 words;
paragraphs median 114 words and 4 sentences. He never writes "note that".

**Captions** open on a noun-phrase headline naming the plotted quantity, terminated by
a period — "Flavor content of the three active mass eigenstates." — never on an
interpretation. Then `{\it Top}:` / `{\it Bottom}:` / `\emph{Left}:`, then the
configuration, then what the curves are, with the visual encoding in parentheses
("(black stars)", "(dashed)"). Close with a pointer. Three to seven sentences.

**Figures**: no gridlines anywhere; ticks inward on all four sides; axis labels as
*quantity, symbol* [unit] with powers of ten folded into the bracket; curves annotated
in place where there is room, legends only for what cannot be; uncertainty as a filled
band; legend boxes with a thin black frame.

**Two standing rules he has stated explicitly.** Never leave an audit undone unless he
says stop or pause. And never use the `", and"` tic in either form — that is criterion
21, which he added.

---

## 3. What was done this session — do not redo any of it

**The paper was written from nothing.** 23 pages, nine figures, four tables, two
listings, five appendices, following the sibling's structure with the method section
replaced entirely.

**Figures.** All nine are produced by notebook 28, which computes Magνs's own numbers
live and reads every other code's from the frozen `notebooks/external_*.json`
datasets. They have been through three rounds of his revisions: house style (tight
axis limits, black legend borders, capitalized labels, plain `1` rather than `10^0`),
then layout (panel counts, orders, one-column widths), then the specific fixes listed
in §5 below.

**Audit criteria already run and cleared:**

- **9**, American English — `travelled` → `traveled`; three other hits were false positives.
- **11 / 19**, short and one-sentence paragraphs — 18 merges.
- **12**, long paragraphs — 7 splits; **0 paragraphs over 170 words remain**.
- **15 / 21**, the `", and"` tic — **44 in the first pass, then 33 more under the
  stricter criterion 21**. Classify by hand; the serial comma and genuine compound
  predicates sharing one subject stay.
- **16b**, voice — measured against the right population, see §2.
- **16e**, captions — all nine rewritten to his pattern.
- **20**, apparatus prose — 11 sentences cut, including the one he quoted at me,
  "Now the part that a paper should not soften."

---

## 4. Measured facts the paper rests on — do not re-derive, but do re-check if you touch them

Every number below was measured this session. Criterion 1 says recompute rather than
re-read; these are given so you know *what* to recompute and what the answer was.

- **The order-two Magnus method is the midpoint slab product.** Not similar — the same
  expression. The one-point Gauss–Legendre node is the slab midpoint. Their probability
  matrices agree to between `1e-15` and `1e-13` over N = 2 to 16 384. `magnus_exp_order`
  1 and 2 both select the same one-node scheme; **the correct name is order 2**. Because
  they are one expression, Fig. 2 draws only one curve for the two — see §7.
- **The tolerance floor is `rtol=1e-12`.** At `1e-14` the residual is an order of
  magnitude *worse* at two and three flavors. At four and five, a `DOP853` reference
  itself moves by `4e-11` between `1e-12` and `1e-13`, which is above what was being
  plotted — hence the per-panel oracle floor now drawn in Fig. 1.
- **The Earth cross-code plane is limited by the Earth models, not the solvers.** Magνs
  floors at `1.6e-7` there while self-converged to `3e-12`; moving the matched potential
  by 1.6 parts in 10⁷ shifts that floor fivefold. The probability inherits the relative
  error of `V_CC` roughly one for one. **The scale is not tuned to minimise the residual**
  — the value used is what the two codes' own constants imply.
- **The adiabatic hybrid answers none of the paper's figures.** `strategy_info` reports
  `engine='cumulative'` for the shock with and without breakpoints, and `engine='magnus'`
  for the tabulated Sun. It does engage on the smooth analytic exponential solar profile
  at loose tolerance. Section 4.8 describes it correctly; no figure exercises it.
- **Standard `P(νμ→νμ)` is flat at 1 above ~100 GeV through the Earth** — measured
  0.9988 to 1.0000 over 1–30 TeV. The vacuum term falls as 1/E. This is why the sterile
  panels sit three decades higher in energy.
- **Term counts** are 1, 1, 2, 3, 5, 9, 17, 33, 65, 129 at orders one to ten. An earlier
  draft had these wrong.
- Packaging: `resources/` reaches neither the sdist (62 entries) nor the wheel (29 files).
  Verified by building both.
- 1207 tests collected. **Coverage of 93% is transcribed from an earlier handover, not
  measured this session** — the one number in the paper I did not verify.

---

## 5. Traps that cost time, so they do not cost it again

- **`hamiltonian_3nu_nsi` returns `V_CC` times the epsilon matrix alone**, and is zero
  when every epsilon is. The standard matter term must be added beside it. Omitting it
  put a reference 0.178 away from the answer instead of 1.5e-5.
- **`|U|²` is indexed final-flavor-first; every `osc_prob` return is initial-first.**
  The transpose is the whole difference and costs exactly the asymmetry of the matrix,
  a few times 1e-5 — which reads as an accuracy difference, not a convention one.
- **`distance_traveled_inside_earth` returns kilometres**; every `osc_prob` baseline is
  in eV⁻¹. Passing the raw value returns a converged, unitary, wrong answer.
- **A figure drawn one column wide still stretches if the float is `figure*`.** Two
  figures were wrong this way for a round.
- **A 3+1 Earth oscillogram at GeV energies costs an hour and aliases anyway** — 21 s
  per row against 0.56 s at three flavors. It is drawn at TeV energies for that reason.
- **`pgrep -f <pattern>` matches the waiting shell's own command line.** Two polling
  loops span forever because of it. Use the harness's background-task notification.

---

## 6. To do

**Run the audit to a clean report.** These criteria have not been swept, or not
finished:

1. **1** — recompute every number. §4 lists what was measured and how; anything not
   there has been inherited. The 93% coverage figure is the known gap.
2. **3, 5, 6, 7, 8** — clarity, omissions, leftovers, internal consistency, and
   consistency against every other section. Counts stated in prose must match what the
   figures draw; several figures changed panel counts late.
3. **14** — never put the pointer before the name.
4. **17** — no self-weakening statements.
5. **18** — every figure and table introduced by a paragraph that *opens* on it. This
   was done for the nine figures; **the four tables were not checked**.
6. **16a, 16c, 16d, 16f, 16g, 16h, 16i** — paragraph shape, orphaned equations, number
   formatting, apparatus, signposting, section titles, habits.
7. **13, dangling words — run LAST and iterate to clean.** It was run once, early,
   against text that has since changed substantially. Detect from justification, not
   paragraph structure; never filter on word count; fix by *adding* about 50−n
   characters, never by cutting. A working detector is:

   ```
   pdftotext -layout main.pdf out.txt
   ```

   then flag any prose line ending in sentence punctuation, shorter than 0.55 of the
   previous line's length, sitting under a full-measure line. Exclude the bibliography.

**Then report**: what was found, what was changed, and what you judged not worth the
prose it would cost. Editorial calls belong to the author — recommend, do not decide.

---

## 7. Figure 2, rebuilt last — the state the paper is actually in

Everything in §§1–6 predates this. Nothing else moved.

**What was wrong.** The right panel measured the truncation orders against a `DOP853`
solution at `rtol=1e-13`. That reference is itself uncertain at `2e-13`, so orders four
and six flattened at about `1e-13` and the panel showed the *reference's* floor while
appearing to show the *method's*. It also drew a separate slab-product curve lying on
top of order two, which is the same expression twice.

**What it is now.** The reference is a midpoint slab product carried at **fifty decimal
digits** in `mpmath` at N = 4096, 8192, 16 384 and 32 768, Richardson-extrapolated
**three times** — the midpoint product's error runs in even powers of the slab width, so
`(4P(2N)-P(N))/3`, `(16Q(2N)-Q(N))/15` and `(64R(2N)-R(N))/63` remove `N^-2`, `N^-4` and
`N^-6` in turn. **Its own convergence is `3.0e-20`**, against a target of `1e-18`.

Two things had to be right for that to work, and only one of them is obvious:

- The Richardson combination must stay in `mpmath`. An earlier version returned `float`
  probabilities from each slab product and combined them in double precision; it
  converged to `1.1e-16`, which is double-precision epsilon and not the extrapolation.
  `mp_reference` therefore returns an object array of `mpf`, and the deviation of each
  Magnus curve from it is taken in `mpmath` too.
- The curves are still double-precision, so they floor where the *composition's*
  round-off puts them. That floor is real and is now visible: order six bottoms at
  `2.1e-15` near 1024 slabs and **rises afterwards**, order four bottoms at `4.5e-14`
  near 8192. Do not read the upturn as an error.

**The measured panel**, which the prose in §4.5 now quotes and which you should not
re-derive:

| N | order 2 | order 4 | order 6 |
|---|---|---|---|
| 128 | 2.83e-04 | 3.98e-07 | 6.11e-10 |
| 1024 | 4.39e-06 | 9.82e-11 | 2.18e-15 |
| 16384 | 1.71e-08 | 1.17e-13 | 1.17e-13 |

Successive ratios are 4.0, 15.9 and 63.7 — `N^-2`, `N^-4`, `N^-6`, as they should be.
Order two never reaches `1e-8` on the panel; order four is below it by 512 slabs and
order six by 128.

**The reference is cached, and the cache is the part most likely to confuse you.**
`notebooks/mpmath_phase_reference.json`, about 3 kB, tracked. Two and a half minutes of
`mpmath` is not worth paying on every rebuild. It is keyed on a SHA-256 **fingerprint of
the configuration** — nine samples of the Hamiltonian along the trajectory, the
baseline, the working precision, the slab counts — so:

- change a mixing angle, the energy, the potential or the baseline and the notebook
  recomputes and rewrites the file by itself;
- change none of them and it prints `reference read from mpmath_phase_reference.json,
  unchanged configuration ac1b2fc890ef` and costs nothing.

**Never invalidate it by hand, and never edit the numbers in it.** If you suspect it,
delete it and let the notebook rebuild it.

**Files touched in this last pass**, all already done:

| File | What changed |
|---|---|
| `resources/paper/main.tex` | Fig. 2 caption rewritten; §4.5 paragraph requoted to the table above |
| `resources/paper/figs/phase_vs_profile.pdf` | Regenerated |
| `resources/paper/README.md` | A paragraph on the cache, under "The figures" |
| `notebooks/make_notebooks.py` | Notebook 28's Fig. 2 cells; `mpmath`, `json`, `hashlib` imports |
| `notebooks/mpmath_phase_reference.json` | New, tracked |
| `tests/test_file_tree.py` | `TREE` entry for that file |
| `pyproject.toml` | `mpmath` added to the `notebooks` optional-dependency group |

**One presentation detail worth not undoing.** The three curves are labelled in place
rather than by a legend, and the labels are rotated to the angle their curve makes *on
the page*, computed from `ax.transData` after `tight_layout` — not from the exponent. A
fixed rotation cannot serve three slopes that differ by a factor of three on a
non-square panel; it put the order-two label through its own curve and off the right
edge. `label_along` in the notebook does this.

**What the audit still has to check here**, since this text is newer than every sweep in
§3: the Fig. 2 caption and the §4.5 paragraph have not been through criteria 20 and 21,
nor through 13. Treat them as unaudited prose.

**Two files in this folder are deliberately untracked**: this brief and
`audit-criteria.md`. They are working documents for the writing, not part of the
package, and `tests/test_file_tree.py` would demand `TREE` entries for them if they were
tracked. `git status` will show them as untracked for as long as they exist. **Do not
"clean them up", and do not `git add` them.**

**The state the gates are in**, checked at the end of this session and not worth
repeating unless you change something they can see:

- `tests/test_file_tree.py` — 7 passed, after `python3 tests/test_file_tree.py --write`
  picked up the new JSON in `docs/source/installation.rst`.
- `ruff` on `notebooks/make_notebooks.py` and `tests/test_file_tree.py` — clean.
- The paper compiles at **23 pages, zero undefined references, one overfull hbox**.
- Notebook 28 executes in **170 s** and all 28 notebooks match the generator.
- The full test suite was **not** run; the author asked for it to be skipped, and
  nothing changed under `src/`.
