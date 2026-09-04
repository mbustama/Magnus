# Pending paper edits

Running list of items raised but not yet applied to `main.tex`. Deliberately untracked,
like `audit-report.md` and `HANDOVER-audit.md` (a tracked file needs a `TREE` entry in
`tests/test_file_tree.py` and a regenerated `docs/source/installation.rst`).

## 1. `N` -> `N_{\rm slabs}` where appropriate

Sec. 4.5 introduces the slab count as bare `$N$`; Sec. 4.3 now writes `$N_{\rm slabs}^{-2}$`.
Pick one. `N_{\rm slabs}` echoes the code's `n_slabs` argument and disambiguates against the
other two things `N` currently denotes.

Slab count -- rename these:

| line | usage |
|------|-------|
| 424  | `N^{-2}` (order-two convergence) |
| 459  | `$N$` slabs (the definition, Sec. 4.5 opening) |
| 472  | `N^{-p}` |
| 484  | `N^{-4}`, `N^{-6}` |
| 486  | `N^{-p}`, `$N$`, `N^{-2}` |
| 488  | `$N$` x3, `\sqrt{N}` |
| 517  | `N^{-p}` |
| 1336 | composing `$N$` constant-density operators |

NOT the slab count -- leave, or give their own subscript:

| line | usage |
|------|-------|
| 528  | `N = 108`, `N = 4096` -- stack size for `eigh` timing |
| 712  | cumulative scan -- number of baselines |
| 1364 | "a further factor of `$N$`" -- nuSQuIDS phase samples |

## 2. Appendix rename, `n` -> `k` for the Magnus order

Sec. 4.1 now uses `k` (matching `\equ{magnus}`'s `\sum_{k\geq1}\Omega_k` and line 366).
The appendix at line 1506 onward still has `\Omega_n`, `n-1`, `n-2` in four places.
Mechanical; `n` is the flavour count elsewhere, which is why it moved.

## 3. Sec. 4.7, the 4x4 closed form

"there being no closed form for a $4\times4$ or $5\times5$ Hermitian eigenproblem"
contradicts the Abel--Ruffini passage in Sec. 1, which grants that the quartic is soluble.
Insert "practical", or split the two cases: the quintic has no solution in radicals at all,
the quartic has one not worth evaluating at these sizes. `expmkernels.py:37` has it right.

## 4. Term counts collide

Paper: `\Omega_7` has 17 terms, `\Omega_{10}` has 129 -- right-nested chains in the lower
`\Omega_m`, matching `magnus.py:1046`. But `expansionterms.py:41` says 26 at order 6, 211 at
order 8, 1918 at order 10 -- fully expanded commutator words in `A` alone. Both correct,
different objects. One clause in either place saying which is being counted.

## 5. "the cost panel"

Fig. 2's caption introduces panels positionally (Top left / Top right / ...) then uses
"the cost panel" once, unglossed, at the end. Body says "the lower left panel". Use the
positional name in both.

## 6. Unitarity numbers: maxima or means?

Sec. 4.7's "4e-16 for a single 3x3 and 4e-15 for a stack of 4096" reproduces as a *maximum*
over the stack -- median and mean are flat at 7.5e-16 / 8.9e-16 from N=1 to N=40960, so
nothing accumulates and the growth is extreme-value sampling. Reworded on that basis.
Still to confirm: whether the probability-output figures (3e-12 -> 1.6e-11 across four
decades in the number of points) are also maxima. If they are means, that clause needs
different wording, since a rising mean would be a real effect.

## 7. Two thresholds in the order-two paragraph

The slab count "stops growing below 0.2 GeV" and the resonance leaves the profile below
0.27 GeV. If these are the same feature, name it once.

## 8. Fig. 2, lower left: drop "Ref." from the legend

`notebooks/make_notebooks.py:13809`

    CODES = [('dop853', 'Ref.: Runge-Kutta order 8 (DOP853)', INK, '-', None, None)]
                        ^^^^^^ drop

to `'Runge-Kutta order 8 (DOP853)'`. Not cosmetic: DOP853 is a *competitor* in this panel,
one of the seven timed configurations, while the referee is the mpmath midpoint slab product
Richardson-extrapolated three times. Labelling it "Ref." is what led the Sec. 4.3 draft to
describe the 50-digit reference as "a stand-alone Runge-Kutta order-8 solver, DOP853 as
implemented in scipy" -- which cannot be right, scipy being double precision.

The comment three lines above says the same thing ("the solver they are all measured
against"); worth rewording to "the solver they are all timed against" so competitor and
referee stay distinct.

Requires regenerating Fig. 2. Cheap relative to Fig. 13 but not free -- batch it with any
other figure edits rather than doing it alone.

## 9. Sec. 4.4: the order-six display is cited to the wrong paper

The three-commutator form the code implements --

    C1 = [a1, a2];  C2 = -(1/60)[a1, 2a3 + C1]
    Omega^(6) = a1 + a3/12 + (1/240)[-20a1 - a3 + C1, a2 + C2]

-- is Eq. (251) of the Blanes-Casas-Oteo-Ros review (`Blanes:2008xlr`, arXiv:0810.5488),
with the alphas its Eq. (257). The review attributes it to BCR, BIT 42 (2002) 262, and
states on p. 96 that three commutators is the *minimum* for sixth order.

`Blanes:2000bit` (BIT 40, 2000) gives a *different* order-six scheme, its Eq. (3.10), in
terms of B^(0), B^(1), B^(2) and needing four commutators. Same nodes, different
arrangement. So the display should cite `Blanes:2008xlr`, not `Blanes:2000bit`.

Add to the .bib (used for the minimality claim, and the order-8 scheme):

    @article{Blanes:2002opt,
          author         = "Blanes, S. and Casas, F. and Ros, J.",
          title          = "{High order optimized geometric integrators for linear
                            differential equations}",
          journal        = "BIT Numer. Math.",
          volume         = "42",
          pages          = "262--284",
          year           = "2002"
    }

Also: "the Gauss-Legendre schemes run out at three nodes" is true of this implementation,
not of the method -- BCR 2000 Sec. 4.1(iii) constructs an 8th-order four-node scheme
explicitly (Eqs. 3.21, 3.22, 4.3; ten commutators, six in the 2002 paper). Word it as a
property of Magnus.

## 10. Future direction: order-eight Gauss-Legendre collocation

Concrete and bounded. Four nodes at v1 = (1/2)sqrt((3+2 sqrt(6/5))/7),
v2 = (1/2)sqrt((3-2 sqrt(6/5))/7); B^(0..3) from BCR 2000 Eq. (4.3); Q1..Q7 from its
Eq. (3.21); assemble by Eq. (3.22). A fourth branch in `_magnus_gl`, a `_GL4_NODES`,
and `MAGNUS_EXP_ORDER_MAX_GL = 8`.

Verification is ready-made: the existing Simpson order-8 path is an independent oracle at
the same requested order, the mpmath reference bounds both, and h^8 convergence is a sharp
test. Payoff is the cost driver Sec. 4.3 already measured -- four evaluations of H per slab
against Simpson's sixty-to-five-hundred, which is what makes order ten cost over a minute
per probability against about 12 ms.

## 11. Order-8 collocation: what is verified, and the R^(4) derivation

Independently verified (Fable, reading BIT 40 (2000) at 300 dpi plus its text layer, and
re-implementing from scratch):

* The order-8 scheme of BIT Eqs. (3.21), (3.22), (4.3) as transcribed in item 10 is correct
  line by line. Local-error slope ~9 across six seed/t0 combinations, and the same with exact
  integrals instead of the 4-node quadrature, so order 8 does not depend on the quadrature.
* Q7 = -(1/42)[B0,[B0, Q3 - (1/3)Q4 + h Q5]]. BIT PRINTS -1/42; there is no typo in the
  coefficient. (An earlier note here claimed otherwise -- that was a misread of the stacked
  fraction at low resolution.) The one real defect in (3.21) is a stray comma, "Q3 - (1/3)Q4,
  +hQ5", which must be a plus.
* (3.21)-(3.22) reproduces the raw R-form (3.11)-(3.20) to 1.4e-16.
* Review Eq. (251) and BIT Eq. (3.10) are both order 6 but differ by O(h^7) on the same three
  nodes -- different constructions, confirming item 9's citation fix.

Alpha notation for order 8 is a DERIVATION, not a transcription: R^(4) appears in neither
paper. Derived exactly (Gauss-Jordan in rationals on T^(4)_ij = (1-(-1)^(i+j))/((i+j)2^(i+j))):

    R^(4) = [[  9/4,    0,  -15,     0],
             [    0,   75,    0,  -420],
             [  -15,    0,  180,     0],
             [    0, -420,    0,  2800]]

so alpha1 = (9/4)A0 - 15 A2, alpha2 = 75 A1 - 420 A3, alpha3 = -15 A0 + 180 A2,
alpha4 = -420 A1 + 2800 A3, with A^(i) = h B^(i). Its even sub-block reproduces R^(3), as it
must. NOT YET VERIFIED: the alpha-form Omega^(8) built on this. Rewriting (3.21) through
R^(4) needs its own convergence check before anything is printed.

Two traps for any hand conversion between the two papers:
* review A^(i) = h * B^(i) (BIT 3.4) -- normalizations differ by one power of h;
* BIT measures node offsets from the slab MIDPOINT, the review from the slab START.

---

## TO DO: re-run the Earth column of Fig. 11 (opened 2026-09-05)

`notebooks/prem_chord_common.py`'s `vcc` evaluated the matter potential with a Python
list comprehension over every position: 6.1 us a point, against 0.063 us for the
vectorized form now in place. Same numbers, to 3.3e-16 -- only the timings differed.

**The handicap fell on one code only.** Magnus is driven through that helper
(`gen_prem_benchmarks.py` passes `prof['vcc']`); NuOscProbExact's Route A
(`append_npe_rtol_prem.py`) builds its potential from `earth.earth_slabs` and never
calls it. Share of Magnus's stored Earth time that was the helper, at 2 flavours:
38% at rtol 1e-3, 44% at 1e-4, 49% at 1e-6, 61% at 1e-8.

**What needs re-running**, on an idle machine:
* Magnus orders 4, 6 and 8, all four flavour cases, Earth chord only
  (`gen_prem_benchmarks.py`, after deleting the three Magnus series per case).
* Then notebook 28, then the figure.

**What does NOT:** NuOscProbExact (never touched the helper); the whole exponential
column (its `vcc` was already vectorized); `prem_chord_reference.json` (the helper
returned identical values, so the reference is sound).

**Paper numbers that move with it**, all in the two paragraphs after "The ceiling is
the flavor count" in `\subsection{A smooth profile...}`:
* "5.9 us per Magnus slab at order 4" and the "0.045 us" it is set against;
* the factor of 130, and the "366 overtakes 130" sentence;
* the ~600 us fixed refinement cost;
* the break-even costs per sample: 6.5 us, 1.7 us, 144 ns.
Unaffected: the 366x slab ratio and the sample counts (972 vs 34,782; 2606 vs 278,494),
which are geometry rather than timing.

Provisional arithmetic, not a re-run: corrected, Magnus at 1e-6 is about 515 us against
NuOscProbExact's 603 us, so it is already ahead there and the crossover sits nearer 1e-6
than 1e-8. Re-measure before printing any of it.
