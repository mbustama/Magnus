# Report: which NuOscProbExact notebooks to adapt, and what to add

**Written:** 2026-08-08. **Status: proposal. Nothing has been written yet — this is the report
the handover asked for, and it waits on your decision.**

Base verified: branch `notebooks`, `a02699b`, clean tree, constants `1.25 65536 True`.

---

## 0. The answer, up front

Six of the twenty are worth adapting. Ranked by what they buy per minute of CI:

| Rank | Adapt | Est. runtime | Why |
|---|---|---|---|
| 1 | 13 Antineutrinos | ~40 s | The convention is written down but never *demonstrated*; verified below that both half-right variants return plausible wrong numbers |
| 2 | 15 Numerical edge cases | ~60 s | Cheap, and the natural home for the nine warning classes |
| 3 | 11 Exact vs approximations | ~120 s | Magνs's version is genuinely stronger than theirs |
| 4 | 12 Ordering and octant | ~150 s | Straightforward, physically useful, both NuFIT sets ship |
| 5 | 09 Performance | ~180 s | Overlaps theirs more than the handover implies — see §3 |
| 6 | 08 Unusual density profiles | ~200 s | Only if the `t_breakpoints` asymmetry is the point |

Plus **four Magνs-only notebooks** (§2), of which I would build **one** first: *"What `rtol`/`atol`
actually promise."*

**My recommendation: do 1, 2, 3 and the `rtol`/`atol` one.** That is ~+6 min on a 37 min budget
and covers the two things a user is most likely to get wrong (the antineutrino convention, and
reading a tolerance as an error bound). The rest can wait for evidence anyone wants them.

Two things I found while checking the proposals are **not** notebook questions and want a
separate decision — see §5. One of them is a silent wrong answer in a public wrapper.

---

## 1. The twenty, re-checked

Source: `/home/mbustamante/Research/NuOscProb/NuOscProbExact/notebooks/` — note the handover's
§2 path is missing the `Research/` component. Structure read from the executed `.ipynb`, not from
their generator, so the cell counts are what actually shipped.

| # | Theirs | Cells / figs | Magνs | Verdict |
|---|---|---|---|---|
| 01 | Basics | 22 / 0 | 01 | covered |
| 02 | Vacuum | 10 / 3 | 02, 03 | covered |
| 03 | Matter, NSI, LIV | 13 / 3 | 02, 03, 08, 09 | covered |
| 04 | Oscillograms | 10 / 2 | 06 | covered |
| 05 | Bi-probability | 15 / 4 | 05 | covered |
| 06 | Earth, PREM, chords | 23 / 4 | partial | skip — 02 §6, 04 and `methodology.rst` cover it |
| 07 | Earth probabilities | 18 / 3 | 04, 06 | covered |
| 08 | Unusual profiles | 15 / 4 | inside 02, 03 | **adapt (6)** |
| 09 | Performance | 24 / 1 | — | **adapt (5)** |
| 10 | Paper figures | 32 / 14 | n/a | blocked — no Magνs paper. Your call |
| 11 | Exact vs approximations | 13 / 2 | — | **adapt (3)** |
| 12 | Ordering and octant | 15 / 4 | — | **adapt (4)** |
| 13 | Antineutrinos | 13 / 3 | — | **adapt (1)** |
| 14 | Solar and MSW | 18 / 2 | 12, 13 | covered, and Magνs's is better |
| 15 | Edge cases | 16 / 1 | — | **adapt (2)** |
| 16 | Four neutrinos | 33 / 2 | 07 | covered |
| 17 | Cross-checks | 18 / 2 | — | **feasible, contra the handover** — see §3 |
| 18 | Evolution operator | 23 / 1 | 11 | covered |
| 19 | Animations | 15 / 4 | — | skip — cost |
| 20 | Arbitrary Hamiltonian | 43 / 6 | — | fold into the Magνs-only "bring your own H" |

### 1.1 Antineutrinos — verified, not assumed

The handover says this is the most defect-prone convention in the package. I built the
demonstration and ran it. At 2 GeV, 1300 km, ρ = 2.848 g cm⁻³, NuFIT 6.0 SK-NO, channel μ→e:

| Construction | P(μ→e) |
|---|---|
| neutrino | 0.074897 |
| **antineutrino, both flips (correct)** | **0.014452** |
| conjugate PMNS only | 0.056774 |
| flip potential only | 0.023478 |

The manual path reproduces `osc_prob_3nu_matter_constant_density(..., nubar=...)` to all six
digits for both signs, so the comparison is against the package rather than against my algebra.
All four numbers are unremarkable probabilities — nothing about the two wrong ones looks wrong.
That is the notebook.

*Correction to the handover:* it says Magνs has "no antineutrino notebook", which is true, but
antineutrinos are not absent — notebook 01 shows `nubar=True` in one cell, and notebook 05
(biprobability) uses `nubar` throughout. The gap is the **demonstration of the convention**, not
coverage of the flag.

### 1.2 Edge cases — the warning classes are the hook

The handover's §4 lists five warning classes. There are **nine**:

`DensityUnitWarning`, `MagnusHighOrderCostWarning`, `ScalarHamiltonianWarning`,
`MagnusConvergenceWarning`, `ToleranceNotAchievedWarning`, `HybridCertificationWarning`,
`UnmarkedDiscontinuityWarning`, `HiddenFeatureWarning`, `PhaseAveragingWarning`.

The last three of the `ToleranceNotAchieved` family and `HiddenFeatureWarning` are missing from
the handover's list. A notebook that shows what each one means, and which are advisory rather
than diagnostic, is worth more than the degeneracy cases their 15 spends most of its length on.

---

## 2. Magνs-only proposals

### 2.1 "What `rtol`/`atol` actually promise" — build this one

The single most valuable notebook I can propose, and it is cheap. A tolerance is a **stopping
criterion**, not an error bound. I ran a case that makes the point without any contrivance —
solar exponential profile, 2ν, 10 MeV, half a solar radius, `rtol=1e-6, atol=1e-8`:

```
convergence_info: n_slabs=20000  n_slab_edges=20000  tolerance_achieved=False
MagnusConvergenceWarning + ToleranceNotAchievedWarning both fire
```

so the notebook has a live, non-toy example of the ladder hitting its ceiling and *saying so* —
next to the ordinary case where it converges quietly. Add `n_slab_edges` vs `n_slabs` (a nominal
2→3 slab step is a 16→17 edge step on an Earth chord) and the PR #35 defect where the ladder
certified an agreement between two nearly identical grids.

### 2.2 "Which engine answered, and why"

`cross_check_strategies` works and is genuinely interesting. The documented example runs four
engines across three families and reports a 2.9e-04 independent spread between `ip_exp` and
`magnus`, with `hybrid` certified. That is a good notebook. **But see §5.2 before building it** —
the function has a failure mode that a notebook would either have to teach or trip over.

### 2.3 "Bring your own Hamiltonian"

Their 20 is the counterpart and is their most substantial notebook (43 cells). Magνs's version
should cover `osc_prob_earth`/`osc_prob_sun` with a user `H_func`, the vectorisation trick
(measured **4.6×**, bit-identical output, per `index.rst`), `ScalarHamiltonianWarning`, and the
palindrome declaration the Earth entry points make for the caller.

### 2.4 "When averaging rescues you and when it does not"

Notebooks 13 and 14 are already the two halves. A short notebook putting them side by side —
phase error falls 53× under averaging, envelope error does not move — would be mostly assembly.
Cheapest of the four; least new information.

---

## 3. Where I read the source differently from the handover

**17 Cross-checks is feasible, not blocked.** The handover dismisses it as "would need an
external code installed". Their notebook 17 has two checks, and only one needs an external code:
the nuSQuIDS comparison reads a **stored JSON**, and the second check codes the
**Zaglauer–Schwarzer** matter eigenvalues inline as a closed form. A Magνs adaptation of the
second half needs no new dependency at all, and would be a real independent oracle — which,
given `crosscheck-cannot-reach-shared-blindness`, is worth more to this package than to theirs.
I still rank it below the six above on cost, but it should not be written off.

**09 Performance overlaps more than stated.** The handover implies the palindrome is something
"Magνs has to show that they do not". Their 09 already has a *"What the palindrome is worth"*
section. Magνs's genuinely distinct material is the batched separable engine and
`BATCH_WORKING_ENTRIES` (1.19–1.38×). The measured palindrome numbers are still worth showing —
0.91× on plain PREM, 1.4–1.67× on an expensive `H_func` — but as a *sharper* version of their
section, not as new ground.

**Minor:** `pytest --timeout` is not available in this environment (no `pytest-timeout`); the
handover's §0 does not use it, but I did and it failed fast.

---

## 4. Cost

Budget is the constraint. Current total is ~37 min; `notebooks.yml` is paths-filtered but any
notebook change runs the whole set.

My recommended four add roughly **6 minutes**. The runtime estimates in §0 are estimates —
anchored to measured comparables in the existing set (08 → 25 s, 05 → 12 s, 04 → 64 s, 11 →
147 s) — not measurements, and I would not defend them to better than a factor of two until the
notebooks exist.

Nothing I propose needs the frozen-oracle pattern. If 17 is ever built, its nuSQuIDS half would.

---

## 5. Two findings that are not notebook questions

Both surfaced while verifying proposals. Neither is urgent; both want your decision separately.

### 5.1 A g cm⁻³ density passed without the unit flag silently returns the vacuum answer

`density_matter_is_in_g_per_cm3` defaults to **False**. Passing a density in g cm⁻³ without
setting it is therefore read as already-in-natural-units, and the matter potential collapses to
nothing:

```
vacuum                            P(mu->e) = 0.046125
rho=2.848, flag omitted           P        = 0.046125   <- exactly the vacuum answer
rho=2.848, flag set               P        = 0.074897
```

No warning fires. `DensityUnitWarning` guards only the **opposite** direction — a density too
*large* to be g cm⁻³, i.e. one converted twice. The under-conversion direction is unguarded, and
it is the one a new user hits, because 2.848 is exactly what a table gives you. It returns a
plausible probability that happens to be the vacuum one.

This is the same shape as the twice-applied antineutrino sign: a self-consistent wrong answer.
I have not proposed a fix; the asymmetry may be deliberate.

### 5.2 `cross_check_strategies` reports `max_spread = 0.0` when it runs nothing

Called with an entry point that has no `strategy` parameter, every engine declines and the
result is:

```
ran = ()      certified = {}      max_spread = 0.0      max_spread_independent = 0.0
```

with no warning and no exception. `osc_prob` itself is such an entry point — and it is the
function notebooks 01, 05, 06 and 07 call directly, so it is the one a reader is most likely to
try. A caller who checks `max_spread` and sees `0.0` reads it as perfect agreement; it means no
engine ran. `ran = [lab for lab in wanted if lab in answers]` is not followed by any empty-set
guard.

My own first call to this function made exactly that mistake, which is how I found it.

Given that this function exists *specifically* to catch methods that certify themselves while
wrong, a silent all-declined path seems worth closing — a warning when `ran` is empty would do
it. That is a src change and outside the notebook brief, so I have not made it.

---

## 6. What I need from you

1. **Which of the six adaptations**, and **which of the four Magνs-only** — or just my
   recommended four (13, 15, 11, and `rtol`/`atol`)?
2. **Their 10, the paper figures.** Is there a paper or a figure you want reproduced? The
   handover flags this as yours to answer.
3. **§5.2** — should I fix the empty-`ran` case, or leave it and teach it in the notebook?
4. **§5.1** — deliberate, or worth a warning?
