# Scoping: the sterile matter projector's single scalar r (issue #47)

Measured at commit `273eb99` on branch `scope-sterile-projector`, 2026-09-06.
Probe scripts live in the session scratchpad
(`/tmp/claude-1000/-home-mbustamante-Research-magnus/a0d2e3ee-05d1-4f5e-8936-3df18ac27060/scratchpad/probe*.py`);
the construction they use is summarized in the appendix so the numbers can be
re-derived without them.

## The short version

1. **The issue understates the error by an order of magnitude.** The 2.1e-2 in
   the warning is the *off-resonance* scale. Near the sterile matter resonance
   the isoscalar default is wrong by up to **0.41 in probability** on a
   core-crossing chord (cos θ_z = −1.0, 2.4 TeV, ν, e-row), 0.35 at the
   issue's own cos θ_z = −0.95, 0.23 in P(ν̄_μ→ν̄_μ) itself. The number is
   converged (flat from rtol 1e-6 to 1e-8 and under n_slabs floors of 64 and
   128) and it is not a spike: the resonance dip shifts ~10% in energy, so
   errors of 0.1–0.36 persist across a ~2-TeV-wide band.
2. **No choice of scalar fixes it.** The best possible scalar at the worst
   point — found by scanning r, not knowable a priori — still leaves 7.4e-3,
   seven times the default tolerance. Band-wide, the best fixed scalar leaves
   ~7e-2. At 3+2 one scalar serves two resonance bands and does worse.
3. **This holds at every Δm²₄₁.** The worst-case error is invariant as Δm²₄₁
   runs over 0.3–6 eV²; only the band moves (E_res ≈ 2.5 TeV × Δm²₄₁/eV²).
   Small mixing does not rescue it: s₁₄ = s₂₄ = 0.05 still gives 8.8e-2.
4. **Off the core, a scalar is fine.** Mantle-only chords (cos θ_z ≳ −0.84)
   with a well-chosen scalar sit at 2e-4 or better. The exposure is
   core-crossing chords in the resonance band, plus a warning hole below.
5. **The full fix is cheaper than the issue implies.** The Earth wrappers
   already break slabs at the PREM layer crossings, Y_e is piecewise constant
   on a subset of those crossings, so a per-layer P is *exactly* representable
   on the existing grid — options (a) and (b) of the issue are numerically the
   same thing on Earth. Runtime cost of P(l) through the general ladder:
   parity (measured 0.85–0.91x, inside this machine's 12–20% drift). The one
   engine that needs real surgery is `separable` — which is the engine that
   answers the common batched-energy Earth call.

The decision reduces to: **the error is large, physical, silent in a
documented-but-leaky way, and unfixable by any scalar on core chords; the fix
is structural but narrower than feared.** Details below.

## Is the issue text current?

Yes. At `273eb99` the projector still takes one `float`
(`matter.py:319`: `proj[k][k] = 0.5*float(ratio_number_neutrons_to_protons)`),
the warning has exactly one call site (`oscprob.py:1741`, reached only from
`_earth_composition`), the Earth wrappers still forward the caller's scalar to
the projector (`oscprob.py:11043-11045`). Nothing recorded in the issue has
been fixed since.

Two corrections to its numbers:

- The measured 2.1e-2 is reproduced at its scale — I find 2.08e-2 in
  P(ν̄_μ→ν̄_μ) at 562 GeV — but it is nowhere near the worst case (see map).
- The warning's own text ("worth about 2e-02 ... twenty times the default
  tolerance") understates the resonance band by ~20x and should be reworded
  regardless of which option is chosen.

## How it was measured

Both arms of every comparison run through the **same** low-level solver
(`oscprob.osc_prob`, the general Magnus ladder) with the **same** density
(the package's own layered-Y_e `rho_func` from `oscprob._earth_composition`)
and the same PREM `t_breakpoints`, so the projector is the only difference:

- **Truth (C):** H(l) = (1/E)·H_vac + V_CC(l)·P(l), with
  r(l) = (1−Y_e)/Y_e from the same `electron_fraction_func_prem` the density
  uses. On PREM this is simultaneously the issue's option (a) and option (b).
- **Scalar arms:** constant P built from r = 1.0 (A, the default),
  r = path-averaged (B, what the warning recommends), r = core's 1.1478 (D),
  r = V_CC-weighted along the chord (E, my candidate better scalar).

Controls, all passed:

- With s₁₄ = s₂₄ = 0 the truth arm matches the scalar arm to **4.4e-16** — the
  P(l) construction changes nothing but the sterile entry.
- The scalar arm matches the shipped `osc_prob_4nu_earth` to 1.4e-6 (wrapper
  at its default tolerance).
- r(l) hits the layer values exactly (crust 1.0194, mantle 1.0173, core
  1.1478 at 10 / 3000 / 6052 km along the cos θ_z = −0.95 chord).
- The headline difference is flat under rtol 1e-6 → 1e-8 and under n_slabs
  floors of 64 and 128 (3.480342e-01 vs 3.480343e-01).
- Engine identity asserted via `strategy_info` / `_engine_probe`: scalar-energy
  Earth calls answer with `magnus` (general ladder), batched-energy Earth calls
  with `separable`, the Sun with `hybrid`. Both arms of every A/B/C/D/E
  comparison use the ladder by construction.

Physics point: package defaults for the 3ν sector (resolved once through
`values_to_unspecified_osc_params` and passed explicitly to every arm),
s₁₄ = 0.15, s₂₄ = 0.10, Δm²₄₁ = 1 eV² unless stated. Full probability matrix;
errors quoted as max over the active 3×3 block unless a channel is named.

## The error, mapped

Worst case over E ∈ [200 GeV, 7.5 TeV] × {ν, ν̄}, per chord:

| cos θ_z | r_pavg | max err, r=1.0 | r=path-avg | r=V_CC-weighted |
|---|---|---|---|---|
| −1.00 | 1.0884 | **4.1e-1** | 1.6e-1 | 8.3e-2 |
| −0.95 | 1.0786 | **3.5e-1** | 1.3e-1 | 6.8e-2 |
| −0.90 | 1.0648 | 2.1e-1 | 8.5e-2 | 4.5e-2 |
| −0.85 | 1.0392 | 1.1e-1 | 4.6e-2 | 2.3e-2 |
| −0.80 | 1.0171 | 2.1e-2 | 2.4e-4 | 5.4e-5 |
| −0.60 | 1.0169 | 9.8e-3 | 2.3e-4 | 4.2e-5 |
| −0.40 | 1.0167 | 4.0e-3 | 1.4e-4 | 3.9e-5 |
| −0.20 | 1.0149 | 1.3e-3 | 1.8e-4 | 5.6e-5 |

Every worst case lands in the resonance band near 2.4–3.2 TeV. The structure
is physical: for ν the error concentrates in ν_e→ν_e/ν_e→ν_s (the s₁₄
resonance), for ν̄ in ν̄_μ→ν̄_μ/ν̄_μ→ν̄_s (the s₂₄ side); τ rows stay at 1e-7.
The r value sets the −V_NC entry, which sets the resonance energy — with
r = 1.0 the ν_e survival dip bottoms near 2.9 TeV, the truth near 3.2 TeV, so
the whole feature is displaced and the pointwise error is order one across the
band.

Scaling:

- **Δm²₄₁ ∈ {0.3, 1, 3, 6} eV²:** worst case 3.57e-1 at every value; only the
  band moves, E_res ≈ 2.5 TeV × Δm²₄₁/eV². Not a corner of parameter space.
- **Mixing:** 8.8e-2 at s₁₄ = s₂₄ = 0.05; 2.6e-1 at 0.10/0.10; 3.5e-1 at
  0.15/0.10; 3.7e-1 at 0.30/0.20.
- **3+2** (s₁₄ = s₂₄ = s₁₅ = 0.10, s₂₅ = 0.05, Δm²₄₁ = 1, Δm²₅₁ = 3 eV²):
  two resonance bands, worst 2.5e-1 for the default, 4.1e-2 for the
  V_CC-weighted scalar. Same disease, more of it.
- **Below ~1 TeV × Δm²₄₁/eV²** (off resonance): default costs 1e-2 to 4e-2;
  a good scalar brings it to 1e-4–2.5e-3.

### A hole in the shipped warning

At cos θ_z = −0.80 the path-averaged r is 1.0171, within the warning's 2%
threshold of the default 1.0, so **no warning fires** — verified. The measured
worst-case error there is **2.1e-2** (3.16 TeV, ν): as large as the figure the
warning quotes for core chords, on a chord the warning declares fine. Near a
resonance the map from composition mismatch to probability error is not
linear, so a threshold on r cannot bound the error. If the warning is kept,
its threshold cannot honestly stay at 2%; if it is tightened to cover this, it
will fire on essentially every up-going 3+1 chord — which is itself an
argument that the warning is the wrong tool.

## Can a better scalar substitute for the fix?

Partially, and measurably not in the place that matters:

- **V_CC-weighting beats path-averaging** about 2x everywhere (it weights the
  core by the potential the neutrino actually feels): worst case 8.3e-2 vs
  1.6e-1, off-resonance ≤ 2.5e-3, mantle-only chords ≤ 2e-4. If the warning
  survives in any form, the value it recommends should be V_CC-weighted rather
  than path-averaged — a one-line change to the warning site.
- **The oracle scalar fails the tolerance test.** Scanning r at the worst
  point, the minimum is 7.4e-3 at r = 1.130 — still 7x the default tolerance,
  found only by knowing the answer. The best r varies with chord, energy,
  channel; there is no formula for it because no single medium produces this
  chord.
- So the mitigation ceiling is: **good everywhere except core chords in the
  resonance band, where nothing scalar gets below ~1e-2.** For 3+1 Earth work
  at the TeV scale — the region these parameters exist to serve — the scalar
  is not an approximation but a different physics problem.

## What would actually have to change

The factorization H_matt(l) = V_CC(l) × P is load-bearing in exactly one
structural place: `osc_prob_matter_std_potential` (and `_nsi`, `_liv`) hand
`(h_vac_energy_indep, VCC_func, h_matt_proj)` **unfused** to the three engine
dispatchers (`oscprob.py:7211/7223/7231`); only the last-resort ladder gets an
opaque closure. Engine by engine:

| Engine | Verdict under P(l) |
|---|---|
| `magnus` (general ladder) | Safe — opaque `H(l)`, samples at quadrature nodes. Answers all scalar-energy Earth calls today. **Measured runtime parity** for P(l). |
| `hybrid` (adiabatic) | Safe — rebuilds its own closure (`oscprob.py:5312-5320`) and hands only `H_of_l` down. Answers the Sun. Declines Earth chords anyway (breakpoints). |
| `cumulative` | Safe — takes `H_func` with energy bound. |
| `separable` | **The real work.** Hoists `mA = -1j*h_matt` once and builds `Vmat = V[:,:,None,None]*mA`, shared across energies (`oscprob.py:4279, 4300-4314`). Answers the batched-energy Earth call — the common production shape. But `Vmat` is already full-size `(n_slabs, m, d, d)`, and the true split is H_E(E) + M(l) with M(l) = V_CC(l)·P(l) energy-independent, so the cross-energy reuse — the actual optimization — survives: build `Vmat` from samples of M(l) instead of scalar-times-constant. Per-layer P is even cheaper: slab edges always include the Y_e boundaries, so P is one matrix per slab (`mA` gains a leading slab axis, the broadcast is unchanged). The physics-informed slab seed `I_V*h_matt` (`oscprob.py:4269`) needs the same generalization. |
| `ip_exp` | Exploits constancy deeply (rotates `h_matt` into the H_E eigenbasis once) but is **restricted to 2 flavors** (`oscprob.py:5037`) — no sterile block can reach it. Ignore. |
| `constant` | Folds `float(VCC_func)*h_matt` (`oscprob.py:4626-4630`); must decline when P is callable. Reachable only through direct `osc_prob_matter_std_potential` calls with constant density plus callable r — a combination the Earth path never makes. One guard. |
| `average` | Fires only when V_CC is constant; same guard as `constant`. |

Also touched:

- `matter.matter_potential_projector` — accept a callable of position (or gain
  a position-resolved sibling); today it `float()`s the argument.
- The three closure sites (std / nsi 4ν+5ν branches / liv) — build
  `htot(enu, l)` with P(l).
- `_PositionProfileCache` (`oscprob.py:2534`) memoises V_CC samples keyed on
  the grid's bytes; P(l) samples want the identical treatment (a sibling cache
  or one fused M(l) cache — the key generalizes, the return shape does not).
- `symmetric_over` mirror-halving on whole-chord requests survives **iff** P(l)
  is built from radius, which the layered Y_e is.
- The `_samples_identical` GL shortcut keeps firing inside layers for a
  per-layer-constant P; a smooth P(l) (the Sun) forgoes it, as the smooth
  density already does.
- Earth wrappers: default the projector to the same `ye_of_r` already closed
  over in `_earth_composition`; retire the warning there.
- Tests/notebooks per the issue's acceptance criteria (3ν bitwise unchanged is
  structurally guaranteed — the sterile block is empty and r never enters the
  3ν Hamiltonian; the frozen datasets need `electron_fraction=0.5` pinned).

**What losing `separable` would cost if it were bypassed instead of fixed:**
measured 7–12x on a 40-energy batched Earth call (0.02 s vs ~0.2 s). An upper
bound — my fallback loop was per-energy Python — but it says the fix must go
through the separable engine, not around it.

## The API question

The owner's standing rule is that arguments passed and returned do not change.

- Option (a)/(b) as scoped here **keeps every signature intact**: the only
  API-adjacent change is `ratio_number_neutrons_to_protons` additionally
  accepting a callable (a type widening, exactly the widening `rho_func`
  already has), plus the Earth wrappers' *default behavior* changing so the
  projector follows the same Y_e the density does.
- The default change moves 3+1/3+2 Earth numbers by up to 0.41 — but that is
  the fix, not a side effect, and the issue's own acceptance criteria embrace
  it: `electron_fraction=0.5` reproduces the old uniform numbers exactly
  (r = (1−0.5)/0.5 = 1 everywhere), 3ν results are untouched.
- Mitigation-only variants (below) change no signature and no default.

## The Sun

The Earth fix does not close the issue: the Sun path takes the same scalar,
fires **no warning at all** (the warning's one call site is Earth-only —
verified), and its default r = 1.0 sits entirely outside the physical range
(~0.14–0.47). A scalar-vs-scalar bracket (r = 1.0 vs the Sun's own ~0.29) on
`osc_prob_4nu_sun` at 10 MeV moves P_ee by **1.3e-2** — consistent with the
issue's 4.5e-3 on an averaged 5ν observable. I could not build a solar truth
arm: the package exposes no Y_e(r) for the Sun (deliberately — the density
profile is already electron number density), so a per-position solar r needs
an SSM composition profile brought in from outside. The Sun's Y_e varies
continuously, which is why only the callable form (a) serves it; the engine
that answers there (`hybrid`) is already opaque to the change.

## Options

**Do the full fix (issue's option (a), which on Earth is also (b)).**
Removes an order-0.1-to-0.4 silent error class at exactly the parameters 3+1
Earth analyses use. Runtime: parity on the ladder (measured), preserved
vectorization on `separable` (by construction — the reuse axis is energy and
M(l) does not depend on energy). Code: the seven touchpoints above; the
genuinely delicate one is `_osc_prob_scan_separable`, second the closure/
dispatch plumbing, third the cache sibling. No signature changes. Risks: this
package's history says dispatch-surface changes breed silent engine-selection
regressions — the cross-engine equivalence tests and an `electron_fraction=0.5`
bitwise pin are the guard rails, and both already exist as patterns.

**Do a mitigation.** Two defensible moves, both small:
(i) change the warning's recommended value (and its internal target) from
path-averaged to V_CC-weighted r — one function, ~2x better everywhere a
scalar helps at all; (ii) reword the warning to state the measured worst case
(0.4 near resonance on core chords, not "about 2e-02") and either drop or
justify the 2% threshold, which is silent at cos θ_z = −0.80 where the error
reaches the very 2e-2 the warning quotes. This leaves every number unchanged
by default and leaves core-chord resonance physics wrong by 1e-2–4e-1 for any
caller who follows the warning's advice. The Sun stays silent and unphysical.

**Leave it and keep the warning.** Defensible only if 3+1/3+2 Earth work in
the TeV band is out of scope for the package's users. The warning as shipped
has the threshold hole, understates the worst case by ~20x, recommends a
scalar that leaves 0.13–0.16 behind near resonance, and does not exist on the
Sun path. If this option is chosen, at minimum the warning text should stop
naming a number it understates.

## What was not measured, and where I am uncertain

- **NSI/LIV Earth variants**: same closure shape confirmed by reading
  (`oscprob.py:7547-7557, 8001-8004`), error maps not measured. No reason to
  expect a different magnitude — the standard projector is added to their
  matrices unchanged.
- **Solar truth arm**: not buildable from the package alone (no Y_e(r) for the
  Sun); only the scalar bracket above. The 4.5e-3/1.3e-2 numbers are
  order-of-magnitude corroboration, not a map.
- **Batched-`separable` vs ladder equivalence** on this exact call was not
  independently re-verified (my arms all ran through the ladder; the package's
  own cross-engine tests own that invariant). The 1.4e-6 wrapper tie-in was
  measured at scalar energies.
- **The 7–12x separable-loss bound** includes per-call Python overhead; a real
  fallback would be somewhat better, and the number drifts with this machine's
  12–20%. Interleaved medians were used for the parity claim.
- The full test suite and the frozen comparison datasets were not run
  (out of scope for a probe; nothing here modified `src/`).
- Convergence of the *truth* arm rests on the same ladder as everything else;
  a genuinely independent external reference (per the comparison notes in
  `docs/source/comparison.rst`) was not rebuilt for P(l). Given the null
  control at 4.4e-16 and the scalar arm's 1.4e-6 tie to the shipped wrapper,
  the exposure is limited to errors common to both arms, which cancel in every
  difference quoted here.

## Appendix: reproducing the truth arm

With the worktree's `src` first on `PYTHONPATH`:

```python
rho  = oscprob._earth_composition(costhz, None, 1.0, None, None, None, None,
                                  'probe', num_flavors=4)   # layered Y_e
VCC  = matter.vcc_func_from_rho_func(rho, 0.0, 1.0, 0.5, nubar, False, True)
r_of_l = lambda l: earth.neutron_to_proton_ratio_from_electron_fraction(
    earth.electron_fraction_func_prem(
        earth.earth_radial_distance_from_depth(costhz, l/gd.UNIT_KM)))
# H(l) = (1/E) H_vac + VCC(l) diag(1,0,0,0) + VCC(l) r(l) diag(0,0,0,1/2)
P    = oscprob.osc_prob(H_func, 0.0, chord, t_breakpoints=prem_edges,
                        rtol=1e-6, atol=1e-6)
```

The scalar arms replace the last two terms with
`VCC(l) * matter_potential_projector(4, r)`. Everything else — density,
breakpoints, solver, tolerances, physics point — is shared.

## An argument for the fix that this scoping did not make

Added 2026-09-06, after checking what NuOscProbExact does with composition.

**The benchmarks are pinned to a uniform `Y_e = 0.5` partly because of this
defect**, and fixing it would let them move to the more physical setting.

NuOscProbExact is not limited to a uniform electron fraction. Its
`electron_fraction` parameter takes a scalar, an array *or* a callable, and the
package ships `earth.electron_fraction_prem`, a PREM-layered profile. What it
does by default is take `gd.ELECTRON_FRACTION_EARTH_CRUST = 0.5`. So the flat
value in Figure 11 and in `tabl{prem}` is a *choice*, not a limitation of the
other code: both packages could be run on layered PREM composition.

The reason not to, today, is this defect. Magnus resolves `Y_e` per layer for the
**density**, but the sterile block of the projector takes one scalar `r` for the
whole chord. A layered `Y_e` therefore makes the density and the projector
describe different media -- which is exactly the failure measured above, up to
0.41 near the resonance. A uniform 0.5 makes `r = (1 - Y_e)/Y_e = 1.0` exactly,
which is self-consistent for the sterile sector and sidesteps the whole problem.

So the comparison is currently insulated from the defect by a setting chosen for
that purpose, and the paper says as much. **Fixing the projector is what would
let the benchmarks use layered composition on both sides**, which is the setting
a reader would assume was used and the one an Earth measurement actually has.
That is a reason to do the fix beyond the size of the error itself: it removes a
constraint on how the comparison can be run.

Decision taken 2026-09-06: keep NuOscProbExact's default uniform 0.5 for now, and
state the value in Figure 11's caption so the choice is visible rather than
implicit.
