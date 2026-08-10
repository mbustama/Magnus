# Handover: notebook 25 as an arbiter between codes

**Written:** 2026-08-10, at the end of a long session. Read §0 and §1 before touching anything.

---

## SESSION 2 STATE (2026-08-10, later the same day) — READ THIS FIRST

Work paused mid-brief for a laptop hibernation. **Nothing is committed**; everything below is
staged in the working tree. No background job is running.

### Done, executing, and staged

* **§11, the supernova shock** (§3 items 1 and 4). `gen_shock_benchmarks.py` writes
  `external_shock_benchmarks.json`; it executes notebook 14's own cells, so there is one
  definition of the profile. Notebook 25 rebuilt: **59 cells, 0 errors, 40.6 s.**
* **§10, solar, corrected** (§3 item 4b) and its speed/accuracy panel (item 3).
* **§5's referee fixed** — see "three findings" below.
* **README + `docs/source/index.rst`** aligned (item 5); axis limits tightened (item 6).

### Three findings that changed what the notebooks say

1. **`prem_referee` in `make_notebooks.py` carried a hand-written `diag(1, 0, 0, 0)`** — a
   *sixth* copy of the matter projector, still pre-fix, invisible to
   `test_matter_projector_is_defined_in_exactly_one_place` because that test compares
   builders and cannot see a copy written out in a notebook. The tell was Magnus and
   NuOscProbExact agreeing with *each other* to 3.7e-04 while both sat 2.6e-02 from the
   referee. It now reads `matter.matter_potential_projector`. **Re-measured since, and the
   verdict inverts:** Magnus goes from 2.573e-02 to **4.508e-08** against the corrected
   referee, NuOscProbExact sits at 2.6e-04 and converges *non-monotonically* (32 slabs per
   segment is worse than 8). Two disciplines in the write-up: 4.5e-08 is below the referee's
   own 4.1e-07 discretisation error, so the claim is "agrees within the referee's
   uncertainty", not an accuracy figure; and the retracted §7 row "PREM, 3+1 ... by ~10^3" is
   replaced by two rows, because cost and accuracy now point in opposite directions there.
2. **This brief's prescribed unitarity diagnostic cannot detect the nuSQuIDS failure.**
   "Sum the flavour probabilities and assert they are 1" passes at **1e-16 on output
   containing P_ee = 3.09**, because nuSQuIDS evolves SU(3) coefficients whose identity
   component *is* the trace — the sum is conserved by construction however badly the
   traceless part integrates. The check that bites is **every probability in [0, 1]**.
   The real cause was the **solver tolerance**: below `rel_error = 1e-6` this calculation
   does not return probabilities. The clustered energy grid and the wrong track were both
   suspected and both **rejected by measurement** (the track was still wrong and is fixed).
3. **The solar speed/accuracy panel this brief specifies is not buildable as specified.**
   Scoring both codes against the analytic adiabatic limit on one axis needs ~1e-5
   resolution; sampling noise in nuSQuIDS's window mean is ~1.5e-2 after 900 evaluations.
   The panel now plots **uncertainty against cost**, with the two entries labelled as the
   different quantities they are. nuSQuIDS's window mean was verified to converge onto
   Magnus at 1.4σ, 0.7σ, 0.3σ (1, 5, 15 MeV, 300 samples each).

### Traps added to §7's list, learned the hard way here

* **A generator whose stdout is its JSON must silence the notebook cells it executes.**
  Notebook 14 prints "ray = 4729 oscillation lengths..." and that one line made the output
  unparseable while every number in it was correct. Fixed with `redirect_stdout(sys.stderr)`.
* **The energy panel allocates `(n_energies, n_slabs, 3, 3)` complex128 = 144 bytes per
  energy per slab.** Unguarded that is 2.1 GB for Magnus and 6.0 GB for NuOscProbExact, and
  the first run was killed by the kernel 25 minutes in with a zero-byte output file. Both
  are chunked to a 192 MB budget, verified **bit-identical** to the unchunked result, and
  phases now checkpoint atomically to `external_shock_benchmarks.partial.json`.
* **`pgrep -f <script>` matches your own shell's command line.** It reported a dead
  generator as RUNNING for several minutes. Match on the PID.

### Measured answers to three questions this brief asks but does not answer

* **§3b, the solar case for notebook 24: the expansion order is NOT reachable.** The brief
  says to check first, and the answer is that `average=True` on the tabulated solar profile
  returns a **bit-identical** result at `magnus_exp_order` 2, 4 and 6 -- `max|P - P|` is
  exactly `0.000e+00` across 40 energies, and the runtime is 0.70 s in all three cases. The
  adiabatic route never propagates, so it never touches the expansion. That case is a
  *finding*, not a panel, and it should be written as one sentence rather than three
  coincident points.
* **§3 item 3, NSI: it does NOT shrink.** NuOscProbExact has `hamiltonian_{2,3,4}nu_nsi`, so
  this is a genuine cross-code comparison. **The conventions match exactly** -- both build
  `V_CC (diag(1,0,0) + eps)`, NuOscProbExact writing the standard piece in as `1 + eps_ee`
  while Magnus adds it separately, which looks like an off-by-one and is not. Constant
  density, 3 g/cm^3, 1300 km, 2 GeV, `eps_ee = 0.10`, `eps_em = 0.05`: the same `eps` into
  both gives P(numu -> nue) agreeing to **1.7e-16**, and shifting either by one puts them
  3.7e-02 apart. Unlike section 6's V_CC comparison there is no offset to correct.
* **`t_breakpoints`, `n_slabs` and `cumulative` reach every entry point through `**kwargs`
  and are absent from every signature.** They do work -- on the NSI shock path
  `t_breakpoints` moves the answer 2.4e-06 and `n_slabs` 8.8e-07 -- and unknown keywords
  *raise* rather than being swallowed, so there is no silent-failure risk. `cumulative=True`
  changes nothing at 61 baselines because `cumulative='auto'` has already engaged. The gap
  is discoverability, not capability, so the cross-code sections must **state which keywords
  were needed**, or the published result is true but not reproducible from the signature.
  Promoting them to named parameters is the wrong fix: ~10 entry points x 3 declarations is
  the same duplication that hid the matter projector.

### All brief items are now built; what each one landed on

| item | where | result |
|---|---|---|
| §3.1, §3.4 shock vs codes | nb25 §11 | residual belongs to a *different code on each front* |
| §3.2 3+1 | nb25 §10b, §12 | solar departure 0.188; shock, Magnus wins both axes |
| §3.3 NSI | nb25 §10c, §13 | conventions agree to 1.7e-16; Magnus to the referee floor |
| §3.4b solar | nb25 §10 | corrected curve + uncertainty-against-cost panel |
| §3b order | nb24 §4 | 1≡2, 3≡4, 5≡6 exactly; solar bit-identical; default 4 confirmed |
| §3c BSM | nb13 §6, nb14 §6 | same eps: 0.014 averaged vs 0.441 instantaneous |
| §3.5, §3.6 | README, index.rst | reach/generality/observables; tight axes |

**The one cross-notebook inconsistency is deliberate and documented in both places:**
notebook 14 uses `D41 = 1 eV^2` and notebook 25 §12 uses `1e-2`. Notebook 14 compares Magnus
against itself, so no referee is needed. Notebook 25 referees both codes with DOP853, and at
eV scale that means resolving 5.9e6 radians over the ray -- of order a day. The eV-scale case
can be computed but not independently validated.

### Open

1. **The gates.** Full suite (running at the time of writing; `src/` has changed three times
   in ways that touch results, so read any failure as a possibly-stale expectation before
   reading it as a regression). Then
   `make clean html SPHINXOPTS="-n -W --keep-going"` -- an *incremental* build is not the
   docs gate and cannot re-find a warning whose file it has cached.
2. **The two-pile commit split described in §0**, which is still the right shape: the physics
   fixes stand on their own and should land first; the notebook work is a second commit.
3. A background task is logged (not started) to validate and document the `t_breakpoints` /
   `n_slabs` / `cumulative` pass-through keywords on the public entry points. Not part of this
   brief; do not fold it into either pile.

This brief is written to stop you repeating work, and more importantly to stop you repeating
*mistakes*. Almost everything that went wrong in the session behind it went wrong in the same
way: the **measurement setup** was wrong, not the code. Every number below that is quoted was
wrong at least once before it was right.

---

## 0. State of the tree

Branch `notebook-25-perf`, cut from `main` at `17c7dd5`. **Nothing is committed** — the last
commit is `fcd1a82`, and everything since is in the working tree. `src/` has changed twice in
ways that alter physics results, so **run the full suite before committing anything**.

```bash
git -C ~/Research/magnus status --porcelain
ruff check src/ tests/ notebooks/make_notebooks.py
python -m pytest tests/ -q -n auto          # ~9 min; 1050+ expected
cd notebooks && python make_notebooks.py --only 25_
```

Uncommitted work is **two independent piles at different stages of readiness, and they should
be committed separately.** Landing them as one change would present unfinished notebook work as
though it were as settled as the physics.

**Pile 1 — physics fixes. Finished, tested, defensible on their own.** Commit these first; they
stand without any of the notebook work.

| fix | evidence |
|---|---|
| sterile neutral-current term, four sites unified into `matter.matter_potential_projector` | 0.29 → 2.4e-04 on a 3+1 PREM chord; 3 regression tests |
| mixing-angle guard, all four vacuum builders (`hamiltonians/_angles.py`) | silent NaN Hamiltonian → `ValueError` naming the parameter; 6 tests |
| `max_n_tpts_per_slab` ceiling tested the wrong variable | that parameter was never validated at all; 3 tests |
| twelve Earth entry points hard-coded `electron_fraction=0.5` | now parameters; verified **bit-identical** at defaults against a pre-change baseline |

The first of these also **retired a claim in `HANDOVER_OVERHEAD.md`** that was telling every
future session not to investigate 3+1 PREM. See §6.

**Pile 2 — notebook 25. Mid-flight; do not land it as finished.** Sections 9 and 10 execute and
the frozen data is generated, but:

* **§10 currently ships the unfair comparison** — Magnus's *averaged* P_ee against nuSQuIDS's
  *instantaneous* one. The corrected data exists (`external_solar_nusquids.json`); it is not yet
  wired in. See §2.
* The shock, 3+1 and NSI cases are not built at all, and the solar speed/accuracy panel is not
  built. See §3.

**The full suite has not run since the mixing-angle guard went in.** `src/` has changed three
times in ways that touch results — the sterile term, the guard, the Earth parameters — so run it
before committing either pile, and expect to read any failure as a possibly-stale expectation
rather than a regression (that is how five of the six failures earlier in the session resolved).

---

## 1. The lessons, which cost the most to learn

**1. A comparison is a setup before it is a measurement.** Six separate results in this session
were wrong because of the setup rather than the code:

| what was measured | what was actually wrong |
|---|---|
| "Magnus is 6.2e-02 from the reference, flat in rtol" | Magnus at NuFIT 6.1, reference at NuFIT 4.0 |
| "Magnus loses 5-30x on smooth profiles" | control read 1.166 — the ratios were not readable |
| "Magnus is 7.2x faster on a scan" | batched Magnus against a **looped** competitor |
| "Magnus is 49-90x slower on PREM" | converged Magnus against NuOscProbExact at 8 slabs |
| "PREM 3+1 is ~1000x slower and inherent" | the matter Hamiltonian was missing a term |
| "the MSW transition is not visible" | the energy window began above the transition |

Before quoting any ratio: is the **physics point** the same, the **potential convention** the
same, the **observable** the same, and is the competitor being run **the way its authors intend**?

**2. Always carry a control, and interleave it.** The house method is in
`implementation_details.rst`: interleave round-robin, report minima, carry a workload the change
cannot touch. A control that returns 1.00 is the evidence the ratios are readable. **Measure the
control by interleaving too** — two *sequential* calls read machine drift and reported 1.133
where interleaving reported 0.986 on the same machine.

**3. Discard the first timed call.** The first Magnus call of a session pays ~0.7 s to compile
the numba kernel. Charged to a twelve-energy batch that is 58 000 us/probability, and it made
Magnus look like one of the slowest codes on the plot. Autorange like `timeit` as well: a
batched call can finish in a few hundred microseconds, too short to time once.

**4. A referee must be independent of what it referees — including in its Hamiltonian.** A
`scipy` slab product refereeing Magnus is an independent *integrator* but not an independent
*Hamiltonian*. When it was built with the same missing sterile term, it "confirmed" Magnus to
1.6e-07 while both were wrong by 0.29. It only broke open when the parametrisations were
compared against a code that did not share the bug.

**5. Clean convergence is not correctness.** Two referees converged beautifully to the wrong
answer: one because uniform slabs straddled PREM's density discontinuities (the tell was
*non-monotonic* convergence), one because `prem_layer_edges_along_chord` returns only the
*interior* crossings and dropping the two endpoints lost 7 km of a 10 830 km chord — 0.065% of
the path, worth 1.8e-03 in probability.

**6. The same three lines in five files is how a defect survives a max-effort review.** The
matter projector was written out in `hamiltonians4nu`, `hamiltonians5nu`, twice in `oscprob`,
and again in the fuzz suite's own oracle. All five agreed with each other and all five were
wrong. They now come from `matter.matter_potential_projector`, and
`test_matter_projector_is_defined_in_exactly_one_place` catches a sixth copy.

**7. An oracle that calls the code it checks cannot catch that code being wrong.** The fuzz
suite's `H_of` builds the projector itself, deliberately, and that is correct design — it was
just encoding the wrong physics. Fixed by writing the *correct* structure out independently,
not by importing the implementation.

**8. Don't remove a section because it is hard; do remove one that is wrong.** The supernova
section was removed at 6.02e-01 error, which was a setup error on my side, not a result. That
was right. But two earlier sections were removed when they could have been fixed. Prefer fixing.

---

## 2. What is done in notebook 25

### Section 5 (PREM, seven codes) — corrected this session
Reads `notebooks/external_prem_speed_accuracy.json`, copied from NuOscProbExact's
`tests/prem_speed_accuracy.json`. **Both codes batched**; the earlier version looped
NuOscProbExact and inflated the ratio ~5x. Workload is theirs: `costhz = -0.9`,
P(numu -> numu), **NuFIT 4.0 NO**, `sin^2(th14) = sin^2(th24) = 0.1`. Magnus needs
`electron_fraction = 0.5 x 1.0001896490` to match their V_CC.

### Section 9 (exponential) — new
Four single-panel figures, 2/3/4/5 flavours, from `external_profile_benchmarks.json` generated
by `notebooks/gen_profile_benchmarks.py`. All codes timed in one process on one machine.
Verdict: **NuOscProbExact is faster at every accuracy it can reach**, and its error floors at
2.5e-11 and then *rises* (32 768 slabs is worse than 16 384). Magnus continues to 2.9e-13.
Magnus wins on **reach**, not speed. Also a probability-vs-energy panel; the two codes agree to
2.7e-11 over 200 energies.

### Section 10 (solar) — landed, but the nuSQuIDS curve must be AVERAGED

**The version in the tree compares Magnus's averaged P_ee against nuSQuIDS's instantaneous one,
and that is not a fair plot.** It was caught in review and the fix is generated but not yet
wired in. `notebooks/gen_solar_nusquids.py` writes `external_solar_nusquids.json`: each target
energy is evaluated at 21 energies spread over ±5% and averaged, which is what a finite energy
resolution means and is the same device `avgprob.averaged_probabilities_numerically` uses. All
samples for all targets go into **one** `EvolveState` — nuSQuIDS evolves the whole energy vector
together, so batching is free and looping would be dishonest as well as slow.

**`external_solar_nusquids.json` AS GENERATED IS WRONG — DO NOT PLOT IT, AND THE FAULT IS IN
THE CALL, NOT IN nuSQuIDS.** It returns `<P_ee>` of **2.55, 58.9, 3.19** at the low-energy end:
probabilities above 1, so normalisation is lost. An earlier draft of this brief blamed nuSQuIDS
for being a GeV-PeV code. **That is wrong** — nuSQuIDS handles low energies perfectly well, and
the defect is in how it is being driven, how the average is taken, or both. Do not repeat that
excuse; find the bug.

Candidates, roughly in order of suspicion:

* **The track.** `nsq.SunASnu.Track(0.05*gd.SUN_RADIUS*units.km)` assumes the argument is a
  *length in km scaled into nuSQuIDS units*. Check what `SunASnu.Track` actually expects — a
  production radius, a fraction of the solar radius, or a path length — and in which units. A
  wrong track is the easiest way to get nonsense out of a correct solver.
* **The averaging window.** `<P_ee>` is a mean over 21 samples spanning +/-5% in energy. If any
  sample is unphysical the mean inherits it; and if the window is too narrow to cover a whole
  oscillation at the low-energy end, the "average" is a phase sample rather than an average.
  The phase varies as 1/E, so a fixed *fractional* window covers very different numbers of
  oscillations across 0.1-20 MeV. Consider a window sized in phase rather than in percent.
* **The initial state and basis.** `Set_initial_state(state, nsq.Basis.flavor)` with
  `state[:, 0] = 1.0`, and the energy vector is passed **sorted** while the state rows are built
  in the original order. That happens to be harmless here because every row is identical, but it
  is exactly the kind of thing to check rather than assume when the output is unphysical.
* **Interactions.** The solver is constructed with `iinteraction=False`. Confirm that is what is
  wanted for a solar propagation and that nothing else needs enabling.

**The diagnostic to run first:** sum the flavour probabilities at every node and assert they are
1 to within the requested tolerance. Any node off 1 localises the failure to an energy, and that
will say far more than the averaged curve does. Do that *before* regenerating anything --
generating 15 minutes of data and then checking it is the wrong order, which is how this was
missed.

**Cost is set by the solver tolerance, not the energy count** — measured 10.9 s per energy at
1e-6, 1.39 s at 1e-4, 0.54 s at 1e-3. That makes the tolerance nuSQuIDS's dial, which is what
gives solar its **speed/accuracy panel** (still to be built): nuSQuIDS traces a curve over
tolerance, Magnus is a single point at 0.66 s, and both are measured against the analytic
adiabatic limit. State plainly that Magnus's residual against that limit is the non-adiabatic
correction rather than an error, so the referee validates and cannot rank below ~1e-5.

### Section 10 (solar) — as landed
Real BS2005-AGS,OP table. **Magnus 40 averaged energies in 0.664 s; nuSQuIDS 12 instantaneous
energies in 130.6 s** on the same file. nuSQuIDS's output spreads 0.79 between neighbouring
energies — it is sampling a ~13 000-radian phase, which is why the instantaneous probability is
not the observable. Referee is the analytic adiabatic limit, which never propagates. The
residual panel is the **non-adiabatic correction** (1e-7 to 1e-5, growing with energy), not an
error — Magnus is the more correct of the two there, so the referee validates and cannot rank.

---

## 3. What is left, in the order it was asked for

1. **Supernova shock vs other codes.** Removed at 6.02e-01 error. **Base it on notebook 14**,
   which has solved this: `ENERGY = 15 MeV`, the probability along the ray as a function of
   distance travelled, and a DOP853 reference frozen in `shock_reference.json` with a
   fingerprint guard that re-derives the profile. My error was reading the JSON's shape instead
   of notebook 14's configuration — check what `Ls` spans (79 911 to 80 000 km in the file,
   which is *not* the full ray), what `L0` it propagates from, and which channel `P` holds.
   NuOscProbExact needs a second driver here: its batched route returns the probability at the
   *end* of a slab chain, and this case wants it at 61 points *within* one — drive
   `evolution_operator_3nu_slabs` and accumulate.
2. **3+1 examples in solar and shock**, timed against other codes.
3. **NSI examples in solar and shock**, timed against other codes.
4. **Probability-vs-energy panels** for the shock (exponential and solar have theirs).
4b. **Solar speed/accuracy panel**, and swap the solar probability panel's nuSQuIDS curve for
   an averaged one — **after fixing the nuSQuIDS call**, since the generated file is wrong. §2.
4c. **Expansion-order section in notebook 24** — three cases, order-labelled speed/accuracy
   curves, and a recommendation. §3b, including the constraint that `gl` gives only three
   distinct orders.
4d. **BSM sections in notebooks 13 and 14** — NSI and 3+1, each against the standard 3nu curve.
   §3c, including the check that NSI is 2nu/3nu only and that the adiabatic route reaches 4nu.
5. **README and `docs/source`**: align the "what Magnus is good for" section with reach /
   generality / pre-packaged observables. See §5 below.
6. **Tight axis limits everywhere.** No dead margin left or right of any curve. Solar is done;
   the section 9 panels and the section 5 panels are not.

---

## 3b. Requested: an expansion-order section in **notebook 24**

Not notebook 25 — this belongs in the performance notebook. Three cases, each with an
**accuracy-against-speed plot whose points are labelled by `magnus_exp_order`**:

1. three-flavour propagation through the Earth (PREM),
2. the averaged solar oscillation,
3. the supernova shock.

Plus prose recommending **how to choose the order, and stating the default**.

### What the order actually does, and the two constraints that shape the section

**The order is never raised automatically.** `oscprob.py:2895` states it: *"osc_prob runs at a
fixed `magnus_exp_order`; the requested tolerance is reached by refining the number of slabs
(and, for the quadrature methods, the points per slab), never by raising the order."* So a
tolerance is met by slabs; the order sets the *rate* at which slab refinement pays. This is
exactly why `MagnusConvergenceWarning` says raising the order will not help a too-wide slab.

**Constraint 1 — with the default `gl` integrator only 2, 4 and 6 are distinct.** Measured on a
PREM chord: order 1 and 2 give identical answers, and 3 and 4 give identical answers, because
the Gauss-Legendre nodes coincide. So the curve has **three** points, not six. Do not plot 1
and 3 as if they were separate settings; either omit them or say why they land on top of their
neighbour, which is itself worth explaining.

**Constraint 2 — `gl` raises above order 6.** `magnus._validate` rejects orders above 6 for
`'gl'`. `gd.MAGNUS_EXP_ORDER_MAX` is 10, but reaching 7-10 requires `trapezoid` or `simpson`,
which changes the *integrator* as well as the order — so those points are not on the same curve
and must be drawn as a separate series if they are drawn at all. `_validate` also warns that
order 7 costs ~2.7x order 6, rising to ~17x at order 10.

### Per case

* **Earth.** Straightforward: `osc_prob_3nu_earth(..., magnus_exp_order=n)` at fixed `n_slabs`,
  refereed as in section 5. Fix the slab count so the plot shows order and not the ladder.
* **Solar.** **Check first whether the order is reachable at all**: `average=True` on a smooth
  profile takes the *adiabatic* route, which does not propagate. The expansion is used only in
  the non-adiabatic patches, so the order may barely move the answer — if so, that is the
  finding and the panel should say it rather than showing three coincident points.
* **Shock.** Blocked on the shock case working at all; see §3 item 1.

### The recommendation to write

Default is **4**, and it is the right default: order 6 buys accuracy only where the slabs are
already narrow enough for the series to converge, and costs about three times more per slab
(`_magnus_gl`: 35.7 -> 5.4 us at order 4 against 102.9 -> 8.7 us at order 6 when the
per-slab-constancy shortcut fires). Raise the order when the profile is **smooth** and the
target accuracy is tight; add slabs instead when a `MagnusConvergenceWarning` appears, when the
profile has a jump, or when the accumulated phase is large. Order 1 or 2 is right only when the
Hamiltonian is constant, where `osc_prob` already forces order 1 because every term past
Omega_1 vanishes identically.

---

## 3c. Requested: BSM sections in **notebooks 13 and 14**

Two new sections in each, and in **every one of the four the standard three-flavour curve is
drawn alongside** — the BSM effect is only legible as a departure from it, and a BSM curve on
its own says nothing about size.

**Notebook 13 (tabulated solar model, averaged probability):**
1. `<P_ee>` with **NSI**, against the standard 3nu curve.
2. `<P_ee>` for **3+1**, against the standard 3nu curve.

**Notebook 14 (supernova shock):** the same two, on the shock profile.

### What to check before building, because these are not all symmetric

* **NSI is two- and three-flavour only in Magnus.** `osc_prob_matter_nsi` builds its matter
  matrix for `num_flavors` 2 and 3; there is no 4nu or 5nu NSI path. So "NSI at 3+1" is not
  available and should not be attempted — say so if it comes up.
* **Does `average=True` reach the NSI path?** The averaging dispatcher branches on whether the
  profile is smooth and position-dependent; confirm the NSI wrapper reaches the adiabatic route
  rather than falling into the numerical-window branch, because the two compute *different
  quantities* (the L/E -> infinity limit against an average over a finite window) and mixing
  them across panels would be the same class of error as the instantaneous-vs-averaged one in
  section 10.
* **Does the adiabatic route support four flavours?** `avgprob.averaged_probabilities_adiabatic`
  is called with a `d`-dimensional Hamiltonian; verify it handles `d = 4` and that the level
  crossings it finds are the physical ones before trusting a 3+1 averaged curve.
* **3+1 solar is a direct beneficiary of this session's physics fix.** The sterile state feels
  neither current, so once the actives' common `V_NC` is removed it carries `-V_NC`. Before the
  fix that term was missing and a 3+1 solar curve would have been quietly wrong — by 0.29 on a
  PREM chord, and solar has a far larger density range. Whoever builds this should re-derive
  that the curve moves when `ratio_number_neutrons_to_protons` changes, as a check that the
  sterile entry is live: see `matter.matter_potential_projector`.
* **The shock is piecewise with declared jumps**, so pass `t_breakpoints` in notebook 14's
  sections exactly as notebook 14 already does elsewhere. A BSM comparison run without them
  would be measuring straddled slabs rather than the new physics.

### The point of each panel

Not "Magnus can do NSI" — the notebooks already establish that. The point is **where the BSM
effect is large enough to matter against the standard curve, and at what energies**. For solar
that is the MSW region; for the shock it is the resonance the shock sweeps through. State the
size of the departure in the text, not just the shape.

---

## 4. How to drive the other codes

**NuOscProbExact** (`import slabs`, `import earth`):
* `slabs.probabilities_Nnu_slabs(H, widths)` takes `(n_energies, n_slabs, d, d)` sharing one
  set of widths and composes the batch in one pass. **Always batch.** 2/3/4 flavours; there is
  no 5-flavour route.
* `slabs.probabilities_Nnu_profile(h_of, baseline, rtol=, atol=)` refines by tolerance and
  raises rather than silently degrading when it cannot meet one within `n_max`.

**nuSQuIDS** — minimum three flavours, but reaches five.
* `nsq.SunASnu(path)` reads a BS2005 table. Three things must be fixed to feed it Magnus's
  copy: its reader needs **uniform columns** (the file has heading lines), it indexes columns
  **by position in the full twelve-column layout** (a trimmed file makes it spline the wrong
  ones), and its **default model path does not exist** in this install.
* `nsq.VariableDensity(x, density, ye)` for an arbitrary profile — this is the route for the
  shock, and it is untried.
* **Conventions**, both verified rather than assumed: length x `0.999999858674`
  (`gd.CONV_KM_TO_INV_EV / 5.0677307162e9`, checked in *vacuum* where no density enters, 3.1e-07
  -> 2.9e-09), and density x **0.99190**, found by scanning for the minimum residual. That
  cross-checks against NuOscProbExact's own 0.99209238 to a ratio of 1.000194 — exactly the
  independently measured Magnus-to-NuOscProbExact V_CC offset. Two unrelated routes agreeing is
  why it is trustworthy.

**GLoBES, Prob3++, NuFast** cannot enter these cases: Earth- or constant-density specific.

---

## 5. What Magnus is for, on the evidence

Worth stating because the notebook kept implying the wrong thing. NuOscProbExact is a
*closed-form* code: where a closed form exists, an exact algebraic solution beats a truncated
series, and that is arithmetic rather than a defect. Constant density, piecewise-constant PREM
and standard 3-flavour are what closed forms are built for.

Magnus earns its place on three axes:

* **Reach.** Slab composition floors at ~2.5e-11 on a smooth profile — past ~16 000 slabs,
  round-off in composing that many products beats the discretisation gain, and more slabs make
  it *worse*. Magnus reaches 1.3e-13.
* **Generality.** Arbitrary `H(t)`: a custom Hamiltonian, a BSM term nobody has diagonalised,
  an interpolated profile. And five flavours, where NuOscProbExact has no route.
* **Pre-packaged observables.** `average=True` returns the adiabatic solar average without
  propagating. Nothing else here offers it; the others must resolve every oscillation first.

---

## 6. Physics fixes made this session, uncommitted

* **`hamiltonian_4nu_matter` and `hamiltonian_5nu_matter` omitted the sterile neutral-current
  term.** Actives share `V_NC` so it cancels; a sterile state does not, and keeps
  `-V_NC = (r/2) V_CC`. Cost **0.29 in probability** on a 3+1 PREM chord, flat in tolerance so
  no refinement revealed it. Now one definition in `matter.matter_potential_projector`, used by
  four sites. **This retired the brief's "PREM 3+1 is ~1000x slower and inherent" claim** —
  post-fix Magnus is ~3x slower at 5e-05 and *faster* at 3.5e-09. See `HANDOVER_OVERHEAD.md`
  §10.1 and §10.6, both marked retracted.
* **All four vacuum builders returned a NaN Hamiltonian** for a sine outside [-1, 1], silently.
  Easy to trigger: the 5-flavour signature interleaves each angle with its phase, so grouping
  the angles puts a *phase* in a sine slot. Now `hamiltonians/_angles.validate_sines`.
* **`osc_prob`'s ceiling check tested `max_n_slabs` while reporting `max_n_tpts_per_slab`**, so
  the latter was never validated. Fixed, with `_refinement_params_rejected` kept in step.
* **The twelve Earth entry points hard-coded `electron_fraction=0.5`,** so matching another
  code's matter convention was impossible through them. Now parameters, defaults unchanged —
  verified bit-identical at 3, 4 and 5 flavours against a pre-change baseline.

---

## 7. Traps specific to this repository

* **`git checkout -- notebooks/` reverts the generator**, which lives beside its output. Use
  `git checkout -- 'notebooks/*.ipynb'`.
* **Stage before running the suite.** `test_tree_matches_git` compares `TREE` against
  `git ls-files`, so an untracked file is invisible locally and red in CI.
* **`fig/` PDFs are gitignored build artefacts** — only 2 of 41 are tracked. The notebook's
  inline PNGs are the durable record.
* **An incremental Sphinx build is not the docs gate.** It cannot re-find a warning whose file
  it has cached; this branch carried a broken `:func:` role for nine commits. Use
  `make clean html SPHINXOPTS="-n -W --keep-going"`.
* **Editing a `code(r'''...''')` cell by string splice eats the `'''),` terminator** when it
  sits on the same line as the last statement. This broke the build three times.
* **`--only 25_` rebuilds one notebook in ~30 s** against ~30 min for all 26.
