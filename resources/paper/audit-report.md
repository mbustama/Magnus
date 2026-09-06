# Audit report — Mag(nu)s CPC paper

Run over `resources/paper/main.tex` against `audit-criteria.md`, supplemented with the
lenses and the refutation discipline from `~/Research/NuOscProbReview/reports/`.
Deliberately untracked, like the brief and the criteria beside it.

Paper state at the end: **24 pages, zero undefined references, three overfull boxes.**
It was 23 pages and one overfull box before; the extra page and two of the boxes are
Fig. 2, which is now square-panelled at full `\textwidth` and takes a page of its own.

---

## 1. Corrections to my own earlier reports

Two findings I reported during the run did not survive checking. Recording them because
the discipline that caught them is the point.

- **"All five tables float unreferenced" — wrong.** I grepped for `ref{tab:`; the paper
  references tables through its own `\tabl{}` macro, which expands to that but does not
  contain the string. All five were referenced all along.
- **"Eight paragraphs exceed 170 words" — wrong.** My counter was treating LaTeX math
  tokens as words. Measured on prose alone, **zero** paragraphs exceed 170 words, which
  is what the previous session's handover claimed.

## 1b. Criterion 1 — the numbers the paper rests on

The handover flagged coverage as the one figure transcribed rather than measured. Both
claims in §6 are now verified against the repository as it stands:

| claim | measured |
|---|---|
| "The repository carries 1207 tests" | **1207 collected, 1207 passed, 0 failed** (14 m 55 s) |
| "covering 93% of the library" | **93%** — 4886 statements, 231 missed, 1742 branches |

The suite was run in full, which the previous session had skipped. Nothing in it broke
under the changes to `make_notebooks.py`, `ode_reference` or `make_H_func`.

## 1c. Figure type sizes — the finding that came out of the author's question

Measured from `main.pdf`: body text renders at **9.96 bp**, captions at **7.97 bp**. Every
figure carried text below its own caption size, and Fig. 2 and Fig. 9 were the worst — in
both, everything except the 8 pt axis labels rendered under the caption.

Two causes, and the second is the one that matters:

- **The notebook overrides `notebooks/matplotlibrc` wholesale.** The rc file asks for
  `font.size 14, axes.labelsize 25, ticks 23, legend 15`; the notebook's setup cell set
  `8 / 8 / 7 / 6.6`, and individual calls went as low as **4.6 pt**.
- **Every figure was silently rescaled on inclusion.** `COL = 3.45` and `WIDE = 7.05` are
  both smaller than the true `\columnwidth` (3.487 in) and `\textwidth` (7.224 in), and
  `savefig.bbox='tight'` cropped them further before `\includegraphics` stretched them
  back. Scale factors read from the PDF content stream: **1.008, 1.016, 1.019, 1.049,
  1.062, 1.116, 1.129, 1.182, 1.337**. A nominal point size was therefore never the
  rendered size, and the discrepancy differed by up to 34% between figures — which also
  means line widths and marker sizes drifted figure to figure by the same factors.

Fixed by setting `COL` and `WIDE` to the true widths, drawing Fig. 6 at the 0.90 of
`\textwidth` it is included at, including Figs. 1 and 9 at the width they are drawn, and
raising 30 size arguments to a floor of 8 pt. After the fix no *whole label* in any figure
renders below caption size; what remains under it is math sub- and superscripts, which
behave the same way in the body text.

## 1d. Criterion 16 --- style, measured against his own seven papers

Measured against `~/.claude/.../memory/mb-writing-style.md`, taken off arXiv:1711.11043,
1901.10087, 1506.02645, 1606.06290, 1808.02042, 1610.02096, 2204.04237.

**A correction first.** The first pass measured against `NuOscProbExact` (3.6 first person
per 1000 words), because the handover argued that genre matters and a code paper should
sit lower than a physics letter. That was wrong: the seven-paper baseline is 8.6, and it
is the one asked for. The memory file now says so, and says that a handover proposing a
different target must be flagged rather than silently adopted.

| marker, per 1000 words | at start | now | his seven papers |
|---|---|---|---|
| first person, body | 2.8 | 5.2 | 8.6 |
| parentheses, body | 6.2 | 7.7 | 27.7 |
| parentheses, captions | 11.2 | 18.8 | 27.7 |
| em dashes, body | 6.0 | 5.3 | 3.7 |
| sentence median | 22 w | 22 w | 18 w |

The parenthesis and first-person targets were relaxed by the author once the honest
conversions ran out: the rest would have meant manufacturing asides and rewriting
statements about the code into statements about us, which costs precision.

**What reading found, that no count could.** A style audit is not a marker count, and
these came only from reading the prose:

- **A broken sentence in the introduction**, two clauses welded with no connector
  ("...are solved exactly by the companion code `NuOscProbExact` [7,8] *evaluates* those
  solutions..."). It predates this session; `git log -S` was used to check before saying so.
- **Five captions carried no visual encoding at all.** `solar`, `shock`, `oscillogram` and
  `prem_plane` never told the reader which curve was which --- a direct miss against his
  caption pattern, and the reason the caption parenthesis count was low. The fix and the
  number turned out to be the same thing.
- **Apparatus prose**: "so we define it before using it", the writer explaining his own
  ordering to the reader.
- **Ornate constructions**, mostly in prose added this session: "and instructively so",
  "the panel carries the same conclusion without leaving it", "each factor being itself
  formed in double precision".
- **Three sentences over 60 words** split, including one carrying a ", and so is" weld.

**A self-inflicted defect worth recording.** A regex that converted ", which is X." to
"(X)." split the math in `$2.8$` at its internal period, producing `(worth a factor of
$2)` and leaving `8$ to it.` stranded. Caught by a `$`-parity check over every
parenthesis. Regex edits on LaTeX must be math-aware, and every such edit since has been
checked by hand.

## 2. Errors of fact found and fixed

| what | where | was | is |
|---|---|---|---|
| Convergence rate under slab halving | §5.1 | `2^{order+1}` | `2^{order}` — the measured 4.0, 16.0, 63.8 are 2², 2⁴, 2⁶ |
| Module count | §6 | "thirteen modules **plus** the hamiltonians subpackage" | twelve plus it; the table has 13 rows *including* it |
| Fig. 1 caption | caption | "the colored curves sit on it" | true at 2 and 3 flavors, a few times above at 4 and 5 |
| Conclusions scaling | §9 | `Φ^0.26` against `Φ^0.88` | `Φ^0.21` against `Φ^0.92`, both codes now tuned to deliver |
| Fig. 6 separation point | §7.4 | "agree bitwise until 12 318 km" | the ramp begins at 12 313 km; first drawn point that differs is 12 329 |
| Fig. 6 difference | §7.4 | "differ by 0.385 at 13 000 km" | a fragile point value — see §4 below |

## 3. Defects in the measurements themselves

These were found by rebuilding Fig. 2 and are the substantive part of the audit.

- **The centre panel was not a fair comparison.** DOP853 was asked for `rtol=1e-10`
  while Magnus was asked for `1e-8`, and the solver integrated its three columns
  separately where one call does. Both are now tuned until they *deliver* 1e-8 against
  the reference, and the solver gets the whole matrix at once, worth **2.5–2.8×** to it.
  The measured gap is 3× at 6 rad and 56× at 952, not the 100× claimed.
- **The left panel measured against a contaminated reference.** Its double-precision
  `expm` is itself wrong by `2e-10` at the top of the range — the same defect the
  previous session fixed in the right panel and left in the left one.
- **Orders 8 and 10 are quadrature-limited at the package default.** They have no
  Gauss–Legendre scheme, so they run on Simpson; at 100 samples per slab the slope shown
  would be Simpson's, not Magnus's.
- **The slope fit was mislabelling curves.** A window reaching to `1e-2` caught the
  plunge between the unconverged regime and the round-off floor and fitted it as slope,
  labelling order 8 as `N^-12`. Corrected window gives −2.00, −3.99, −6.01, −8.03,
  −10.05, −12.27.

## 4. Fig. 6, the centre panel

The paper's explanation was right and under-supported. The width parameter widens *both*
fronts, and the one below the forward shock is the contact discontinuity at 12 348 km.
Whether widening it matters is set by its width against the **local oscillation length**
there, measured at 16 km.

Sweeping the width and the energy together is the test that could have refuted this, and
it passes: the columns stay together while the ramp is thinner than that length and part
once it is not, and the width at which they part tracks the oscillation length with
energy (6.3 km at 5 MeV, 15.7 at 15, 46.6 at 45).

The "0.385 at 13 000 km" was not reproducible as drawn: the probability oscillates with a
~16 km period on a 17.5 km grid, so the nearest plotted point shows 0.166. Replaced with
statistics that do not depend on where you stand — a mean difference of 0.22 through the
shocked region, mean survival 0.91 against 0.72.

## 5. Language criteria

| criterion | outcome |
|---|---|
| 9, American English | clean |
| 11, one-sentence paragraphs | 23 found; 16 were headings, listings or section-final. One merged |
| 12, long paragraphs | clean (0 over 170 prose words) |
| 14, pointer before the name | clean |
| 15, "And" openers | 8 found, all fixed |
| 16a, paragraph shape | median 86 words against his own 114 — the paper runs short. Reported, not forced |
| 16c, orphaned equations | clean — all 21 displays have a running lead-in |
| 16d, number formatting | one 4-figure mantissa (2.848 g cm⁻³); defensible, left |
| 17, self-weakening | clean |
| 18, floats introduced | five tables now opened on; `tab:prem`'s introduction moved before its float |
| 20, apparatus prose | 11 fixed, in the body, the abstract and three table captions |
| 21, the ", and" tic | 132 classified by hand; 9 split, the rest legitimate |
| 13, dangling words | iterated to clean: **no prose tail under 16 characters** |

## 6. Left undone, deliberately

- **Paragraph length.** Bringing the median from 86 to his 114 would mean merging
  perhaps thirty paragraphs. That is a voice decision, not a defect.
- **Tails of 16–17 characters.** Three remain. The measure is 61 characters, so these are
  a quarter of a line — inside the zone the criterion says to decide about rather than
  clear.
- **Fig. 2 now takes a whole page.** That follows from square panels at full
  `\textwidth`. Reducing to `0.9\textwidth` would give the page back.
- **Criterion 2, quotations.** The paper quotes no sources verbatim; citation accuracy
  against the originals was not re-checked.

---

## PARKED — Fig. 2, to finish on an idle machine

Two items, both measurement rather than layout, deliberately left for a quiet machine
because they are timing-sensitive.

1. **Order two misses the panel's stated accuracy at four of twenty-five energies.**
   At $\Phi = 41$ to $127$ it delivers $1.01$ to $1.21 \cdot 10^{-8}$ against a $10^{-8}$
   target: its $N^{-2}$ convergence runs into the slab ceiling and no setting on the
   ladder delivers. It *does* deliver at the hardest energy, so the failure is not
   monotone in phase --- the error oscillates with the probability. **Decided: raise the
   order-two slab ceiling and re-measure**, rather than caption around it.
2. **More energies, for smoother curves.** The cost panel carries 25; the curves are
   visibly segmented. Raise the count once the ceiling change is in, so both are measured
   in one pass.

Everything else in Fig. 2 is done: the 2x2 arrangement, the density panel, orders eight
and ten, the energy axis with phase above it in powers of ten, the split legend (solver
in the cost panel, the six configurations with their measured rates in the panel that
measures them), and type sizes set against the page rather than the canvas.
