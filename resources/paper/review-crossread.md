# What the oscillation-probability review offers the Mag(nu)s paper

A cross-read of `~/Research/NuOscProbReview/paper/main.tex` (7 798 lines, 366 bib
entries) against `resources/paper/main.tex` (30 pp., 202 bib entries), done in three
passes. Nothing here is implemented. No text is lifted from the review; every
recommendation is phrased as what our paper would need to say in its own words, and
where a claim of the review's is quoted as motivation it is marked as theirs.

## How the comparison was made, and one correction

My first citation comparison matched **bib keys** and reported 350 review references
absent from ours. That is wrong and I discarded it: our bibliography uses different
keys for the same works (`GLoBES` and `Huber:2007ji` where the review has
`Huber:2004ka`; `Prob3pp` where it has `Prob3`; `nuCraft` where it has
`Wallraff:2014qka`). Matching on **normalized title and arXiv eprint** instead gives
74 shared and 292 genuinely absent. Two claims I would otherwise have made are dead
on that recount: we already cite the standard Magnus review (Blanes et al.) and we
already cite GLoBES, Prob3++ and nuCraft. Everything below survived the recount.

---

# Pass 1 — inventory, section by section

| Review section | Our corresponding place | First read |
|---|---|---|
| 7.2.1 Slab products | 4.5 Composition | We already prove order 2 = midpoint slab product. Their layer-placement result is sharper than our version of it |
| 7.2.2 Parametric resonance | 6.3 Earth / Fig. 5 | We show oscillograms and never say what makes their structure |
| 7.2.3 Magnus integrators | 4.1, 2.2 | They cite the 1990 origin and a 2016 speedup paper; we cite neither |
| 7.2.4 Other structure-preserving | 4.4, 4.7 | Cayley/Crank--Nicolson/Wei--Norman context we lack |
| 7.2.5 Adiabatic and Landau--Zener | 4.8 hybrid | Attribution and the solar/supernova line |
| 7.2.6 Exact solvable profiles | 5.2 engine 3, 5.5 | An exact benchmark exists for exactly our exponential engine |
| 7.2.7 Turbulence | 5.6 Limitations | A frequency-domain criterion for our known blind spot |
| 7.3--7.4 Approximating eigenvalues / amplitude | 4.7 | Mostly closed-form territory; one claim in our favor |
| 7.5 Averaging | 4.9, 5.2 engine 1 | We already handle degeneracy; partial averaging we do not |
| 7.6 Reducing to the observable | 5.6, 7 | Where engine accuracy stops mattering |
| 9.2.2 Magnetic resonance | 4.1, Conclusions | Sixty years of technique for choosing the composition |
| 10.2 What the codes do | 2.2, 8 | 18 codes on a common axis; we name 6 |
| 11 Open problems | 7, Conclusions | One open problem our code is built to answer |

---

# Pass 2 — triage

**Killed.** Four candidates did not survive checking against our own text:

- *Degeneracy in the averaged limit.* I was going to recommend it. Sec. 4.9 already
  carries the block form, Eq. (averaged_blocks), reducing to the singleton case.
- *The Blanes Magnus review.* Already cited, different key.
- *GLoBES / Prob3++ / nuCraft absent.* Already cited, different keys.
- *"Truncation is exactly unitary."* Already ours, and better stated than in the review.

**Held but downgraded.** Broadening the code survey (below, I1) is worth doing but is
housekeeping, not substance; it should not displace anything else.

**Survivors** are the twelve below. Each is graded by whether it needs a new
measurement, since a recommendation that needs one is a different proposition from a
recommendation that needs a citation.

---

# Pass 3 — the recommendations

## A. Provenance we are missing

**A1. The Magnus expansion's arrival in neutrino physics (needs: one citation).**
Sec. 2.2 says the expansion has been used for neutrino oscillations before and that we
do not claim the idea. That sentence currently cites nothing at the point where it
most needs to. The line begins with D'Olivo and Oteo (1990), `DOlivo:1990xs`, and runs
through Ioannisian and Smirnov (which we do cite). *Risk: none. This is the one gap I
would fix first — a generosity sentence that names no one reads as a formality.*

**A2. Casas, Ioannisian et al. (2016), `Casas:2016asi` (needs: one citation, and a
sentence of positioning).** "Efficient numerical integration of neutrino oscillations
in matter." It is the closest prior art to what we do and it reports speedups over
general-purpose ODE solvers of up to two orders of magnitude. We currently measure
ourselves against `solve_ivp` and against other neutrino codes without acknowledging
that this comparison has been made before in the Magnus setting. *Risk: it invites the
question of what we add over it — which Sec. 2.2 can answer (arbitrary flavor count,
arbitrary H, the engine machinery, orders to six), but should answer explicitly.*

**A3. Average Hamiltonian theory (needs: citations, plus a paragraph in Conclusions).**
What we call the Magnus expansion is average Hamiltonian theory in magnetic resonance,
arriving there in 1968 (`Waugh:1968`, `Haeberlen:1968zz`, `Evans:1968`) and becoming a
*design discipline* (`Haeberlen:1976`): sequences chosen so that named orders of the
series vanish identically, with symmetry cancelling classes of higher-order terms. The
review's observation is that the expansion reached neutrinos in 1990 and the design
discipline never did. For us the analogy is exact in one respect and inexact in
another, and the paper should say both: a spectroscopist controls the Hamiltonian,
whereas a matter profile is given — what we control is the *composition*, meaning where
the slab edges fall and how each slab is quadratured. We already exploit half of that
(edges on discontinuities). The symmetry half is untouched. *This is the strongest
outlook paragraph available to us and it costs no new work.*

## B. Sharpening a claim we already make

**B1. Layer placement against integrator order (needs: a measurement we can run).**
Sec. 4.3 says a slab straddling a discontinuity "degrades locally to low order," and we
report the population number 7.8e-4 to 1.3e-12 on declaring edges. The review's
measurement is stronger and more useful: with edges placed uniformly, orders two, four
and six all stall at the *same* floor, so raising the order buys nothing whatever until
the discontinuities are respected. Ours is a statement about one order; theirs is about
the ordering of the two effects. We have the machinery to measure our own version —
error against slab count at orders 2, 4, 6, with and without declared edges, on a PREM
chord. *This is the single most useful new figure panel available. Do not cite their
number; measure ours.*

**B2. Pre-asymptotic behavior (needs: reading our existing data).** The review notes
that below roughly twenty layers every scheme but the slab product converges more
slowly than its nominal order. If our own convergence data show the same, one sentence
in Sec. 4.6 would stop a user concluding from a coarse run that a high order is not
delivering. *Cheap. Check before asserting.*

## C. An exact benchmark for the interaction-picture engine

**C1. The exponential profile is exactly solvable (needs: citations; optionally a
validation test).** Sec. 5.2 introduces the interaction-picture engine for a two-flavor
exponential profile, and the package's own documentation says an exact reference exists
for that case — but the paper cites nothing for it. The two-flavor exponential profile
solves in confluent hypergeometric (Whittaker) functions: `Petcov:1987xd`,
`Toshev:1987jw`, and `Petcov:1987zj`, which we already have. The linear profile solves
in parabolic cylinder functions (Zener, which we have); the tanh profile in
hypergeometric ones. The reduction that produces all of them is in `Notzold:1987cq` and
`Kim:1987ss`, and Balantekin showed the solvable family is fixed by a shape-invariance
condition from supersymmetric quantum mechanics, `Balantekin:1988aq`. For a profile in
none of these classes, `Akhmedov:2008nq` gives a computable perturbative series.
*Two payoffs. The engine gets provenance for the case it was built for, and Sec. 5.5
gets a candidate validation against a closed-form special-function reference rather than
against our own converged run — which is a stronger check than anything we currently
have for that engine.*

## D. Making our own figures say something

**D1. Fig. 5, the oscillogram (needs: citations, and two sentences).** Our introducing
paragraph is entirely computational: chord geometry, mandatory slab edges, the electron
fraction convention. It never says what a reader is looking at. Two physical statements
are available and both are standard: much of the structure on core-crossing
trajectories is parametric enhancement, a resonance between the oscillation phase and
the layering rather than an MSW resonance (`Akhmedov:1988kd`, `Krastev:1989ix`,
`Akhmedov:1998ui`, with the Earth treatment in `Akhmedov:2000js` and the oscillogram
literature in `Akhmedov:2006hb`, `Akhmedov:2012ah`); and in our sterile panels the
matter-driven resonance sits near E ~ dm2_41 / 2|V_s|, which locates the feature to
within tens of percent. *Criterion 18 says the introducing paragraph says what the
float is for, and 16e says a caption must interpret rather than list. This paragraph
currently fails both, and it is our most-looked-at figure.*

**D2. Fig. 9, astrophysical flavor composition (needs: citations).** The section cites
production physics only (Margolis, Stecker, Kelner) and nothing of the flavor-ratio
program the figure actually plots: `Learned:1994wg` for the original argument,
`Athar:2000yw` and `Beacom:2003nh` for its development into the measurement program,
`Farzan:2008eg` for how far the decoherence justifying the averaging can be pushed, and
`Song:2020nfh` for the inverse problem of inferring the production mechanism from the
composition at Earth. *A section that plots flavor composition at Earth and cites none
of this will read as unaware of it.*

## E. Scope statements that strengthen rather than weaken

**E1. Where engine accuracy stops mattering (needs: one paragraph; no measurement).**
We sell a reach of 2.9e-13. The review's Sec. 7.6 measures what sits above a probability
engine in an atmospheric analysis: interpolating an oscillogram off a grid finer than
production practice carries ~2.4e-2 RMS, and ~4.9e-3 after a realistic energy response,
against published engine accuracies of 1e-9 to 1e-6. Our paper should state plainly
where its precision is and is not the binding constraint — it is binding for
cross-code validation, for method development, for regimes with no closed form, and for
averaged observables computed once; it is not binding inside a binned atmospheric fit.
*Criterion 17: a scope statement is precision, not concession. Stating this makes the
rest of our accuracy claims more credible, not less.*

**E2. Phase averaging is not energy smearing (needs: one or two sentences).** Our
`average=True` removes oscillating terms at fixed H. Smearing over reconstructed energy
is a different operation because H itself depends on E, so it moves between
Hamiltonians rather than selecting among one Hamiltonian's frequencies. The review
quantifies the difference as a matter effect — negligible in vacuum, appreciable near a
resonance in matter. A user who takes our averaged output for a detector-resolution
average is making an error we can head off in two sentences.

**E3. The baseline is not a single number (needs: one sentence; `Gaisser:1997eu`).**
Our Earth entry points take L from the chord geometry. An atmospheric neutrino is
produced somewhere in a column tens of kilometers deep, which spreads L by about a
factor of fourteen between vertical and horizontal. nuCraft carries that distribution;
we do not. *An honest scope line, and it costs nothing.*

**E4. No derivatives (needs: one sentence; `Fernandez-Menendez:2025qbi`,
`Granger:2026qdr`).** Three of the eighteen codes the review tabulates now return
gradients — analytically in CHIC, by automatic differentiation in MANGO and nuTens.
We do not, and a fitting user will want to know before adopting. Natural future work.

**E5. Partial averaging (needs: one sentence; `Blennow:2005yk`, `Kiers:1995zj`, and
`Ohlsson:2000mj`, which we have).** We provide complete averaging. The intermediate
case — finite coherence, wave-packet separation, finite resolution — is handled by
damping factors at the assembly stage and is outside what we do.

## F. A claim in our favor that we do not currently make

**F1. We never form a characteristic polynomial (needs: a measurement).** The review's
central accuracy problem is that closed-form SU(n) evaluations divide by the derivative
of the characteristic polynomial at each root, and lose digits as two eigenvalues
approach. It measures an evaluation that avoids forming the polynomial staying below
6e-11 and flat in the gap, where closed forms degrade to about 1e-6. Our general ladder
never forms it either, and our exponential backend reduces powers by Cayley--Hamilton,
which is that same family of construction. *This is a structural advantage we have and
never claim. It must be measured on our own code before it is stated — sweep a
near-degeneracy and show our error is flat in the gap. If it holds it belongs in
Sec. 7; if it does not, we learn something more important than the sentence.*

## G. Turning our vaguest limitation into guidance

**G1. Stochastic and turbulent profiles (needs: citations, and a rewrite of one
limitation).** Sec. 5.6 records that broadband roughness is invisible to the structural
detectors, which is true and unhelpful — a user cannot act on it. The supernova
turbulence literature supplies the criterion in the frequency domain: decompose the
fluctuating part of the profile in Fourier modes, and a mode drives transitions when
its wavenumber matches an eigenvalue splitting, with long-wavelength modes suppressing
the transition instead (`Patton:2013dba`, `Patton:2014lza`; the effect was first studied
for random supernova fluctuations by `Loreti:1995ae`, and its consequence for the shock
imprint is in `Fogli:2006xy`). Restated that way, our limitation becomes advice: the
structure that defeats a spatial-feature detector is the structure whose power sits near
a level splitting, and a user with such a profile should set a slab floor from the
splittings rather than trusting the ladder. *This is the recommendation I would rank
second after A1, because it converts a known weakness into a usable instruction.*

## H. A number only our code can produce

**H1. The cost of the standard solar reduction (needs: a measurement we can run).**
Nearly every solar analysis replaces the three-flavor problem by a two-flavor one at a
rescaled electron density. That reduction is exact in vacuum and approximate in matter;
the literature puts the neglected terms at a few times 1e-3 at solar-core density and
10 MeV, with a next-to-leading correction available and generally unused
(`Kuo:1989qe`, `ParticleDataGroup:2024cfk`, `Lim:2002iz`). We compute the full
three-flavor solar problem, so we can measure the cost of that reduction directly across
energy rather than cite it. *A cheap, genuinely new number, in a section we already
have (6.4), that differentiates us from the codes built on the reduction.*

## I. Housekeeping on the code survey

**I1. Breadth (needs: citations).** The review tabulates eighteen public codes on one
axis — the stage at which each stops being exact. Sec. 2.2 names six. Worth adding at
least `OscProb`, `Kallenborn:2019ilo` (CUDAProb3), `Neurthino`, `Fong:2022oim`
(NuProbe), `NuPert`, and the three differentiable ones from E4. The "where does it stop
being exact" framing is also a cleaner organizing idea than our current prose, though
adopting it wholesale would be a larger edit than this paper needs.

**I2. PEANUTS is our solar comparison, not nuSQuIDS (needs: a citation, possibly a
measurement).** `PEANUTS` implements precisely the adiabatic-plus-decoherence solar
treatment that our averaged solar engine computes, and the review recommends it as the
comparison point for a new implementation. Sec. 8 compares our solar case against
nuSQuIDS. PEANUTS is the closer match and the more informative comparison.

## J. An open problem we are built to answer

**J1. How error accumulates along a composed profile (needs: a measurement; optional
for this paper).** The review lists as open, and says it is unsettled anywhere, whether
the per-step numerical penalty of a composed evaluation adds in proportion to the number
of steps, grows like a random walk, or largely cancels — noting that the integrator
literature bounds the *discretization* error, which is a different quantity. We compose
thousands of slabs, we already track unitarity drift, and we have an mpmath reference.
Measuring the composed error against slab count at fixed discretization would answer a
named open problem at modest cost. *I would not hold the paper for it. It is the
strongest candidate for a follow-up, and worth one sentence in the Conclusions naming
it as open even if we do not do it.*

---

# Ranking

If only five are taken: **A1** (the 1990 origin), **G1** (turbulence as guidance),
**D1** (make the oscillogram say something), **B1** (order versus edge placement,
measured), **E1** (where our precision binds).

If a sixth: **C1**, because an exact special-function benchmark for the exponential
engine is a stronger validation than anything that engine currently has.

Three need new measurements before anything is written: **B1**, **F1**, **H1**.
One needs only a check against data we already have: **B2**.
Everything else is citations and prose.

---

# Implementation status (2026-09-01)

All twelve recommendations are in the manuscript. Three needed measurements, which were
run first; the numbers quoted below are ours, not the review's.

- **B1.** Measured on a PREM chord at `cos th_z = -0.85`, 5 GeV, against a 200 000-slab
  order-six reference. With uniform edges, orders two, four and six stall together near
  `1.2e-4` and stop improving with slab count; at 200 slabs order two is the most
  accurate of the three. With the sixteen boundary crossings declared, the same counts
  give `2.8e-8`, `2.1e-13`, `7.1e-14`. Placed in Sec. 4.6.
- **F1.** Random Hermitian Hamiltonians with the two closest eigenvalues driven together,
  against a 60-digit `mpmath` reference: error flat across ten decades of relative gap,
  median `4e-14` at accumulated phase `1e3` and `2e-12` at `1e5`. Placed in Sec. 4.7.
- **H1.** Full three-flavor solar against the standard two-flavor reduction at a density
  rescaled by `cos^2 th_13`, same profile and parameters: the reduction is low by
  `1.8e-4` in the median over 0.1-20 MeV, rising to `1.9e-3` at 20 MeV, median `8.7e-4`
  above 5 MeV. Placed in Sec. 6.4.

**B2 was dropped, not implemented.** It predicted that low slab counts would converge
more slowly than the nominal order. Our own B1 data show the opposite: with the edges
declared, orders two, four and six sit on their predicted slopes from 50 slabs upward,
and order six is already at the reference floor there. The claim did not survive its own
check, so nothing was written.

**One setup error worth recording.** The first H1 run omitted
`density_is_of_number_of_electrons=True`, so the solar profile was read as a mass
density; the survival probability came back flat in energy, which contradicted our own
Fig. 6 and is what caught it. The same trap is already in
`magnus-comparison-setup-before-measurement`.

41 references were added (202 to 243 entries). Four carry no year, matching the
convention the manuscript already uses for software repositories; they have no year
upstream either, so none was invented.
