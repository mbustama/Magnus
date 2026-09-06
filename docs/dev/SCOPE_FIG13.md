# Scoping: the proposed Figure 13 (averaged solar probability, speed against accuracy)

**Verdict up front.** The quantity does not support the figure as proposed. On the
BS2005-AGS,OP profile at 0.1–20 MeV there are **no non-adiabatic crossings at any of the
four flavour counts**, so `average=True` returns two eigendecompositions with
`P^cross = I`: a single point per flavour count, exact to 1.4e-15 against a 60-digit
evaluation of the same closed form, with **no parameter that moves the returned value at
all** — measured, not argued (§1). There is no accuracy axis to sweep and therefore no
Magnus curve to draw in Figure 11's shape. What the quantity does support is §5.

Everything below was measured in this worktree (v1.1.0, commit 3075e19), with the engine
asserted via `strategy_info` / `oscprob._engine_probe` in every arm: the averaged calls
all answered from the `'average'` engine, the coherent ones from the certified
`'hybrid'`. Machine load stayed at 0.35–0.76 throughout; every timing quoted is an
anatomy with order-of-magnitude margins, not a race.

Setup, common to all probes: the notebook 25 §10 configuration exactly — the
`bs05_agsop.dat` table interpolated linearly in log n_e and clamped at both ends, the ray
l = 0 to R = 0.983 R_sun (the table's last node), NuFIT 6.1 NO, d = 2 params
`{sth: s12, Dm2: D21}`, the 3+1 point of notebook 25 §10b
(s14 = s24 = sqrt(0.10), D41 = 1e-5), and for 3+2 additionally
s15 = s25 = sqrt(0.10), D51 = 3e-5 (my choice; no shipped 3+2 solar point exists — see
§6). Direct `avgprob` calls were validated against the shipped wrapper before use:
max |direct − wrapper| = 0.0 exactly, at every flavour count.

## 1. Question 1: no accuracy curve exists

**Windows.** `find_nonadiabatic_windows` at its defaults finds **zero windows** at every
(d, E) probed — d = 2, 3, 4, 5 × E = 0.1, 0.3, 1, 2, 5, 10, 15, 20 MeV. The adiabaticity
parameter never comes close to the 0.1 threshold:

| gamma_max      | 0.1 MeV | 1 MeV  | 5 MeV  | 20 MeV |
|----------------|---------|--------|--------|--------|
| d = 2          | 4.2e-8  | 4.8e-6 | 9.8e-5 | 5.6e-4 |
| d = 3          | 4.1e-8  | 4.7e-6 | 9.8e-5 | 5.6e-4 |
| d = 4          | 9.6e-7  | 8.7e-5 | 6.7e-4 | 2.9e-3 |
| d = 5          | 8.6e-7  | 6.4e-5 | 4.7e-4 | 2.0e-3 |

The margin is a factor of 34 at its narrowest (d = 4, 20 MeV). With the v1.1.0
composition-following neutron-to-proton ratio instead of the scalar default, gamma_max
rises to 8.8e-3 at d = 4, 20 MeV — still a factor of 11 under threshold, still zero
windows. The `undecided` and `undecided_between_crossings` lists are empty everywhere:
every eigenvalue pair is fully decohered over this ray at all four flavour counts
(including the D41 = 1e-5 and D51 = 3e-5 pairs), so the averaged expression is valid and
un-warned.

**Knob sweep.** With zero windows the structure of `averaged_probabilities_adiabatic`
predicts that nothing can move the value, and the sweep confirms it bit-for-bit. At
(d = 3, 5 MeV), (d = 3, 20 MeV) and (d = 4, 20 MeV), each of

- `n_points` 51 / 801 (default 201),
- `threshold` 0.3 / 0.02 (default 0.1),
- `n_probe` 50 / 1600 (default 200),
- `fd_step_frac` 1e-5 / 1e-8 (default 1e-6),
- `magnus_exp_order` 2 / 4 (default 6),
- `integration_method` trapezoid / simpson (default gl)

returns max |ΔP| = **0.0 exactly** against the default call. `n_points` feeds only the
report flags (the phase integrals are computed after P is formed and never touch it);
the other five act only inside windows, of which there are none.

**The one lever that is not a dial.** Forcing `threshold` below gamma_max does open a
window — a single one covering 40–100% of the ray, because the exceedance is a broad
plateau of small gamma, not a localized crossing. The answer then moves by 1e-8 to
2.7e-7, and the patch usually reports `converged = False` (it is being asked to
coherently propagate most of a 13 000-radian trajectory inside its 32 768-slab cap,
which is exactly the case its own docstring says should decline). An uncertified twitch
at the 1e-7 scale is not an accuracy axis. Two further structural facts close the
question: the patch tolerance is an internal constant (`patch_atol = 1e-7` in
`_local_evolution_operator`), unexposed through `level_crossing_matrix`; and the public
wrapper forwards **no parameters whatsoever** to `averaged_probabilities_adiabatic`
(oscprob.py line 3752 calls it with defaults only), so a user of `average=True` could
not reach a dial even if one existed.

**Cost anatomy** (d = 3, 5 MeV, per energy): full call 16.0 ms, of which the window
search that finds nothing is 11.5 ms (72%), the decoherence report 1.9 ms (12%), and the
two eigendecompositions that are the answer 0.021 ms (0.13%). The fraction of the cost
living in the crossings is zero; ~84% goes to certifying there are none. Per-energy
wrapper cost across flavour counts: 9.3 / 14.0 / 24.0 / 39.5 ms at d = 2 / 3 / 4 / 5
(40 energies at d = 3 is 0.56 s — the shipped "about 0.7 s" claim reproduces).

## 2. Side finding: the shipped 1.3e-05 is a readout convention, not a non-adiabatic correction

Notebook 25 §10 (and `comparison.rst` line 238) scores `average=True` against an
"analytic adiabatic limit" that decoheres at production and **reads out in vacuum**, and
attributes the 1e-7 to 1.3e-5 residual to "the non-adiabatic correction … Magνs carries
it, because its adiabatic route uses the exact crossing probabilities". Measured: **no
crossing probabilities were used** (zero windows, §1), and replacing the reference's
vacuum readout with the matter eigenbasis at the clamped endpoint — n_e(R)/n_e(0) =
8.1e-6, the table's last node — reproduces `average=True` **bit-for-bit (0.0)** at every
energy. The residual is entirely the difference between two endpoint conventions:
readout in matter at 0.983 R_sun, where the trajectory ends, versus readout in vacuum,
where a detector sits. Its growth with energy (1.3e-7 at 1 MeV to 1.3e-5 at 20 MeV)
tracks the matter perturbation V/Δ ∝ E at the surface, not any crossing physics.

Two consequences. For this figure: an accuracy axis drawn "against the analytic
adiabatic limit" measures a convention, not an error — it cannot serve as Figure 13's
y-axis. For the shipped material: notebook 25's lower-panel prose ("it is the
non-adiabatic correction", "below about 1e-5 Magνs is the more correct of the two") is
wrong as written and should be corrected separately; the extension of the ray to true
vacuum (or the statement of the endpoint convention in the caption) is the actual choice
being made there.

## 3. Question 2: the reference is trivial, and it was demonstrated

With the window list empty, the returned object is |V(l0)|² · |V(l1)|²ᵀ — two Hermitian
eigendecompositions and a product. In mpmath (`mp.eighe` at dps = 60, float64
Hamiltonian entries converted exactly) this reproduces the d = 3, 5 MeV call to
**1.38e-15**, with row sums exactly 1, in **1.7 ms** — cheaper than the 13 ms float64
package call, because the reference does not hunt for windows. No reimplementation of
`level_crossing_matrix` is needed. A reference builder for this quantity should (a)
assert the window list is empty and record gamma_max against the threshold, so the
regime is certified rather than assumed, and (b) evaluate the closed form at high
precision. If a profile with genuine crossings is ever wanted, the reference would
additionally need a high-precision integration of the evolution across each window (an
mpmath Magnus or RK integrator over the patch — days of work, and pointless on this
profile).

## 4. Question 3: the competitor gate, measured — and the panel it would have joined is gone

NuOscProbExact has no route to the averaged object; the question was whether a
self-computed mean of its coherent output converges to it, and whether the converged
value depends on the averaging window. Both halves are now measured at d = 2, 5 MeV,
full ray:

**It converges — onto the window-smeared observable.** The mean of `average=True` over
a ±w energy window (201-point grid, uniform in 1/E) shifts deterministically from the
point value: +4.3e-5 at ±2%, +2.7e-4 at ±5%, +1.1e-3 at ±10%, +4.4e-3 at ±20% — the
decohered curve is sloped and curved through the MSW transition, so its window average
is not its point value. Sampled means of the *coherent* probability land on that smeared
limit, not on the point value: NuOscProbExact's slab route (Hamiltonian stacks built
from Magnus's own pieces, so no potential-convention gap; self-converged to ~1.5e-6 by
N = 131 072 slabs; unitarity 1e-13) gives, at 6 144 samples per window,
mean − smeared = −0.12 SEM at ±2% and +0.04 SEM at ±10%; a 48-sample control with the
certified hybrid as sampler agrees. So the plan's gate condition 1 passes and gate
condition 2 fails *as framed* — the converged value depends on the window — but for a
physical reason: the estimator computes the smeared observable, and the window is a
statement about detector resolution, not a numerical knob.

**Why it still cannot carry a curve.** Comparing fairly means smearing Magnus
identically (milliseconds, closed form, noiseless), after which the competitor's
distance from the target is pure sampling noise ∝ 1/√N. Measured cost at converged slab
count: 3.3 ms per sample, so SEM 3.4e-3 costs ~20 s per point, 1e-3 costs ~4 min,
1e-4 costs ~7 h **per energy per flavour count** — against Magnus's single point at
1.4e-15 for ~14 ms. The reachable decades do not overlap; there is no shared axis on
which both draw curves. This is the same conclusion notebook 25 §10 reached for
nuSQuIDS ("sampling noise … three orders of magnitude above the scale the axis would
need to resolve"), now with first-party numbers for NuOscProbExact. Drop it from any
solar panel; the honest caption sentence is that no other code computes this object,
and that a sampled mean of a coherent code converges to the *smeared* observable at
1/√N.

## 5. What the quantity does support — the decision

**(a) Build as proposed: no.** There is no Magnus accuracy curve (§1), no meaningful
competitor curve (§4), and the y-axis the plan inherited from notebook 25 measures a
convention (§2). A Figure 11-format column would contain four isolated points and no
curves.

**(b) Build in a different shape: yes, two candidates.**

- **b1 (recommended).** One column, the profile on top, then the *observable itself* —
  ⟨P_ee⟩(E) over 0.1–20 MeV at 2, 3, 3+1, 3+2 — with the cost per curve (0.4–1.6 s for
  40 energies, measured per-energy 9–40 ms) and the reference agreement (1.4e-15 against
  the mpmath closed form) stated in the caption as numbers, not axes. This is the
  four-flavour-count extension of the existing `fig/solar_averaged.pdf` +
  `fig/solar_3plus1.pdf`, in one figure, and it shows the thing the section actually
  claims: an observable the other codes do not offer, at closed-form cost. No accuracy
  panel — a flat line at 1.4e-15 is a caption sentence, not a plot.
- **b2.** Keep a speed–accuracy plane but as *points*: notebook 25's
  `fig/solar_speed_accuracy.pdf` already is this figure (Magnus star, nuSQuIDS circles,
  axes explicitly declared incomparable). If it is promoted into the paper, the star's
  y-value should first be re-scored against the mpmath reference of the object it
  computes (1.4e-15), not against the vacuum-readout limit (1.3e-5) — see §2.

**(c) "The quantity does not support this figure": true for the proposed shape.** The
finding the figure was meant to carry — a tolerance dial traded against runtime — does
not exist for this quantity, in either code. That is not a gap in the measurement; it is
the point: the averaged solar probability is a closed form, and a closed form has no
convergence knob. Saying that in one caption sentence is stronger than a panel of
points pretending to be curves.

## 6. Notes for whoever builds it

- **The 3+2 parameter point needs blessing.** I used s15 = s25 = sqrt(0.10),
  D51 = 3e-5 to complete the scan; nothing shipped defines a solar 3+2 point. Windows
  stay at zero there (gamma_max 2.0e-3 at 20 MeV), so the §1 conclusion is insensitive
  to the choice, but a published curve is not.
- **The sterile-sector ratio treatment moves the 3+1 answer by 2.2e-2 at 5 MeV**
  (1.0e-5 at 20 MeV): scalar `ratio_number_neutrons_to_protons=1.0` (what notebook 25
  §10b plots) versus the v1.1.0 composition-following callable built from the table's
  hydrogen fraction (n/p = (1−X)/(1+X), 0.47 at the centre to 0.14 at the surface).
  Windows remain zero either way, but any 4/5-flavour curve must state which it uses,
  and the physical choice is the callable.
- **The endpoint convention decides the readout basis** (§2): the ray ends at the
  table's last node with the density clamped, so `average=True` reads out in matter at
  0.983 R_sun, and the two conventions differ by up to ~1.3e-5 at 20 MeV (measured).
  Extending the profile to true vacuum would move the readout basis to vacuum and
  should recover the vacuum-readout limit, but that extension was not probed here.
  Captions should state the endpoint; neither choice is an error.
- **Nothing in the public API reaches the averaged path's parameters** — a caption or
  text must not imply a tunable tolerance for `average=True` on a smooth profile.
- Probe scripts (throwaway, session scratchpad, not committed):
  `probe1_windows.py` (windows/gamma/engine/wrapper-equality),
  `probe2_knobs_attrib_ref.py` (knob sweep, §2 attribution, mpmath reference, cost
  anatomy), `probe6_q3_competitor.py` and `probe7_window_dependence.py` (§4). All
  configurations are stated in full above; each rebuilds from `bs05_agsop.dat` and the
  shipped builders in a few dozen lines.
