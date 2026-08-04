# Handover: does any of this reach a real user? Physically-motivated profiles

**Written:** 2026-08-04, at the close of the session that produced the four robustness commits on
`dev-robustness`.

**Where to work: branch `dev-robustness`, which is where everything below already lives.**
Nothing is pushed and there is no PR.

---

## 0. Verify the base before starting

```bash
git -C ~/Research/magnus log --oneline -1        # must be d647c7a
git -C ~/Research/magnus status --porcelain      # must be empty
python -c "import sys; sys.path.insert(0,'src'); import magnus.oscprob as o; print(o.HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS, hasattr(o,'cross_check_strategies'))"
# must print: 8 True
```

If the seam prints 25, or `cross_check_strategies` is missing, you are on `main` or on an older
commit and **nothing below applies**.

Branch state, so you are not surprised:

| | |
|---|---|
| commits ahead of `main` | 4: `ac65a00`, `3397fcc`, `5a66f01`, `d647c7a` |
| tests | **751 passed**, ~17 min |
| `ruff check src/ tests/` | clean. `ruff check .` reports 63, all pre-existing in `notebooks/` |
| docs | `make html SPHINXOPTS="-W --keep-going"` succeeds |
| bit-identity | 1 of 11 workloads moved vs the start of the work, justified in `FINDINGS_ROBUSTNESS_PROGRAMME.md` §12.1 |

**Read first:** `FINDINGS_ROBUSTNESS_PROGRAMME.md` §§11–12, then this. §12.8 is the list of what
is still open; this brief is one item from it.

---

## 1. The job, and why it is the right next thing

Three tranches of robustness work produced findings ranging from *wrong by 0.54 in probability*
to *four orders of accuracy left on the table*. Asked whether a real user would ever hit them,
the honest answer was: **unknown, and the evidence leans towards no.**

Every catastrophic number came from one of two places:

* a **deliberate adversarial construction** — an unmarked density step, a Gaussian narrower than
  the probe grid, a sawtooth;
* a **random Fourier sum** — the fuzzers' smooth-profile generator, which is a mathematically
  convenient way to make a smooth positive function and has no physics in it at all.

Against the profiles the package actually ships for, the record is clean:

| population | silent misses |
|---|---|
| 164 Earth + solar configurations, d = 2…5, standard/NSI/LIV, ν and ν̄ | **0** |
| 42 workloads on solar, multi-resonance, noisy and a resolvable bump | **0** |
| the same 42, on random Fourier sums | **2** |

So the whole robustness case currently rests on profiles nobody would compute. **The job is to
build a population that a referee would accept as physically motivated, and re-run the existing
instruments against it.**

**A negative result is a real result here.** If the physical families show no silent miss, that
is the answer to the reachability question and should be written up as such — not treated as a
failed hunt, and not "fixed" by making the profiles more adversarial until something breaks.

---

## 2. The profiles to build

Ordered by how much I expect them to pay. Each needs a one-line justification in the code saying
*where the physics comes from*, because a finding on a badly-built profile is a finding about
your construction.

### 2.1 Tabulated density with interpolation kinks — build this first

The most common way a real user's profile stops being smooth, and the least glamorous. A user
loads a density table (a stellar model, a coarse PREM, a simulation snapshot) and interpolates.
Linear interpolation gives a **C⁰-but-not-C¹ kink at every node**; a cubic spline is C² but rings.

Construct: take the package's own solar exponential, sample it at N nodes, and interpolate
linearly. Sweep N so the kink spacing crosses the probe grid: N = 20, 50, 200, 1000, 5000
against `n_probe0 = 200` and `max_n_probe = 6400`.

Why it should pay: the previous session measured a kink (C⁰ not C¹) at **1.448e-02, silently**,
before the fixes. This asks whether ordinary interpolation reproduces that.

### 2.2 Supernova shock front

The textbook physical case where a real density profile is near-discontinuous, and a well-studied
MSW problem (forward shock, reverse shock, contact discontinuity).

Construct: base profile ρ(r) ∝ r^−2.4 (or r^−3), with a forward shock at r_s where the density
jumps by a factor of ~5–10 over a width w, and a reverse shock behind it. Sweep w from 1e-2 down
to 1e-6 of the trajectory. **Two narrow features on one trajectory** is the case the window-merge
logic and the hidden-feature scan should both be exercised on.

Cite a reference in the docstring (Schirato & Fuller, or Fogli et al. on shock effects in SN
neutrino oscillations) so the shape is defensible rather than invented.

### 2.3 Density fluctuations with a Kolmogorov spectrum

Also physically motivated for SN envelopes: δρ/ρ drawn from a power-law spectrum. Unlike a random
Fourier sum with a handful of modes, this deliberately has **power at every scale, including
below the probe grid** — which is exactly the regime `find_hidden_features` was built for and the
one place I expect it to fire on something real.

### 2.4 Earth with a non-PREM crust

`osc_prob_earth` passes PREM's layer edges as `t_breakpoints`, which is why Earth has been clean
throughout. A user with their own crust model, a 3-D tomographic slice, or an added sediment
layer may well *not* pass breakpoints. Construct PREM plus two or three extra crustal layers and
run it **without** `t_breakpoints`.

### 2.5 A real tabulated solar model

If you can get one (BS05/AGSS09-style), use it directly rather than the analytic exponential.
This is the highest-credibility profile available and costs nothing but the download.

---

## 3. What to run — reuse, do not rebuild

Everything you need exists under `docs/dev/adversarial_batteries/` and every script is listed in
its `README.md`. **Each of these takes a profile list; swap yours in rather than writing a new
harness.**

| script | what it answers, once you give it your profiles |
|---|---|
| `fallback_quality.py` | every applicable engine scored on the same request. This is what moved the seam from 25 to 8 |
| `warn_fp.py` | true/false-positive rate of every warning, split-oracle so it is affordable |
| `resolution_fp.py` | the resolution test's FP rate, swept over sub-intervals |
| `crosscheck_acceptance.py` | whether a cross-check between engines would see any defect found |
| `constants_audit2.py` | re-run the constants against the physical population if anything moves |
| `bitident.py` | run before and after any change; 11 workloads that must not move without a reason |

`harness.py` already encodes every API trap (see §5) — build profiles with its helpers.

---

## 4. Pass criteria — write these down before running anything

1. **(P1) No silent miss.** On every physical family, `strategy='auto'` at the default tolerance
   is inside 1e-3 **or** warns. A silent miss is the headline result and must be attributed:
   which engine answered, which mechanism, and whether `t_breakpoints` cures it.
2. **(P2) No new false positives.** `find_hidden_features` must stay at 0 false positives on the
   *smooth* physical families (it is 0 of 67 on the synthetic ones, at every sampling density).
   A physical profile that trips it without a real sub-grid feature is a defect in the detector.
3. **(P3) Detection rate on the families that do hide something.** For the shock and turbulence
   families, the fraction of genuinely-hidden features the scan reports **is the number that
   answers the reachability question.** It is 68–90 % on synthetic Gaussians; if it collapses on
   physical shapes, that matters more than any of the synthetic numbers.
4. **(P4) The seam change holds.** `HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS` moved 25 → 8 on a
   cost measurement over **three** profiles. Re-measure the cost on the physical population; if
   the cumulative scan is not cheaper there, the change needs revisiting.

---

## 5. Traps — every one of these cost real time in *this* session

**Process discipline** (I broke all of these at least once)

* **`pgrep -f foo` matches your own shell's command line.** `pgrep -f "fallback|bitident" | kill`
  killed the command issuing it. **Kill by PID from a targeted `ps`.**
* **Never spawn an `until … sleep … done` waiter.** I left one polling a file for **3 h 38 m**
  after killing the job it was waiting for, and at one point had three polling the same file.
  **Poll directly and sparsely.**
* **One heavy job at a time.** Three concurrent measurements under `ulimit -v 10000000` on a
  15 GB box got two of them OOM-killed with no message — the failure mode is a kill, not a
  slowdown. Check `free -g` before starting a second.
* **One output file per run.** Two jobs writing the same path produced a stale read that looked
  exactly like a real test failure, and I chased it.
* Foreground commands are killed at **2 minutes**; use `run_in_background` for anything longer.
* Always `(ulimit -v 10000000; …)`.

**Measurement discipline**

* **Absolute times under load are worthless.** Alternate and carry a control the change cannot
  touch; if it does not return ~1.00×, discard the round. Mine returned 0.99× and that is the
  only reason the seam numbers are usable.
* **A stalled battery is usually the oracle.** `solve_ivp` cost explodes at low energy and high
  flavour count. `warn_fp.py` died at 100/180 the first time for exactly this; it now splits the
  oracle — `expm` (exact, free) for piecewise profiles, `solve_ivp` only for smooth ones at
  ≥ 30 MeV. **Do the same for your profiles.**
* **A population that does not contain the workload you are about to change is not evidence
  about it.** This bit me on `threshold0`: a fixed-baseline sweep said a rule was safe, and an
  energy scan said it was 20× worse. Include points, baseline scans *and* energy scans.

**Physics / API**

* `earth.density_matter_func_prem(r)` takes a **radius from Earth's centre in km**; its second
  argument is `tol`, not `costhz`. Route through
  `earth.earth_radial_distance_from_depth(costhz, l/gd.UNIT_KM)`.
* `matter.vcc_func_from_rho_func`'s **7th positional argument is
  `density_is_of_number_of_electrons`, not `nubar`.** Use `harness.vcc_of`, where both are
  keywords and the mistake is unrepresentable.
* A `rho_func` must return a **scalar for scalar input**; `harness.scalarize` handles it.
* A generator that drives a profile negative gets `ValueError: rho_func must be non-negative` —
  that is **correct validation**, not a defect. Guard the loop and count them.

**Docs**

* `sphinx -W` treats any new warning as an error, including inside a `jupyter-execute` block.
* An RST **simple table's column width must fit its widest cell** — a 10-character heading in a
  9-wide column failed the docs build this session.

---

## 6. Do NOT rebuild these — they were built, measured and rejected

Full reasoning in `FINDINGS_ROBUSTNESS_PROGRAMME.md` §§11.2, 11.3, 12.2.

| idea | why it is not there |
|---|---|
| `strategy='auto'` cross-checking itself in the weak band | Fired **0 times** on 200 profiles. With `GAMMA_TO_ERROR` mis-calibrated 2× it *still* fired 0 times, because the trigger was computed from that same constant — self-referential. Verifying every window-free result also fires 0 times. What is left there is engines wrong **together** |
| `threshold0` derived from the tolerance | Single point 13711× better, **energy scan 20× worse** |
| `MagnusConvergenceWarning` keyed to the returned refinement level | Firings 70 → 53 but **true positives 17 → 4**. The mechanism is kept private as `magnus._deferred_slab_norm` |
| auto-inserting `t_breakpoints` at a detected feature | Partial cure (3–46×) and it changes dispatch. Choosing a grid is the caller's call |

---

## 7. Where I expect this to be harder than it looks

1. **Building a defensible profile is itself the work.** A supernova shock has a shape, a speed
   and a width, and getting them wrong produces a finding about your construction. Budget for
   reading before coding, and cite the source in the docstring.
2. **The oracle will fight you.** SN profiles over long baselines at low energy are exactly where
   `solve_ivp` is unaffordable. Keep spans short, energies ≥ 30 MeV, and use `expm` wherever the
   profile is piecewise-constant — it is *exact* there, not an approximation.
3. **"Physically motivated" is a spectrum.** Be explicit per profile about where it sits: a real
   tabulated model is not the same evidence as an analytic shock with plausible parameters. Say
   which is which in the write-up.
4. **The interpolation-kink case (§2.1) is the most likely to pay and the least interesting to
   build.** Do it first anyway.
5. **Expect the answer to be "mostly no".** The Earth/solar surface has been clean across 164 +
   42 configurations. If the physical population is also clean, the correct conclusion is that
   the package's exposures need a user-supplied non-smooth profile — which is worth stating
   plainly, and changes how the warnings should be pitched.

---

## 8. Deliberately not done in this session

* **The 42-workload `fallback_quality.py` re-run after the seam change.** Killed partway to save
  time; the two silent misses were verified directly instead (1.19e-03 → 7.26e-08 and
  1.68e-03 → 1.87e-09) and the mechanism is deterministic. **Re-run it before merge.**
* **`ENGINE_FAMILIES` groups the cumulative scan with the general ladder.** Right about shared
  blind spots, wrong about accuracy — they differ by four orders on the same request. Documented
  as a known limit rather than changed, because splitting them would let two engines that *do*
  share a blind spot vouch for each other.
* **Nothing is pushed and there is no PR.**
* **GitHub Pages is still disabled**, so the "Documentation Deployment" workflow fails on every
  commit. Settings → Pages → source "GitHub Actions". Nothing in the codebase can fix it.
