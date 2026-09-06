# Plan: a new Figure 13 -- averaged probability on the BS2005-AGS,OP model

**Not started.** One column, in Figure 11's format: speed against accuracy for the
*averaged* solar probability, at each flavour count, on the tabulated
BS2005-AGS,OP profile.

## What the plotted quantity is, and why that decides the figure

`average=True` does not return a time-average of an oscillating probability. On a
position-dependent profile it calls
`avgprob.averaged_probabilities_adiabatic`, which

> decoheres in the eigenbasis **at production**, transports along the levels of
> the instantaneous Hamiltonian (with the exact crossing probabilities wherever
> the evolution stops being adiabatic), and reads out in the eigenbasis at
> detection.

So it is a closed-form object costing an eigendecomposition at each end plus
crossing probabilities -- not an integration, and not a sampled mean. That is why
it runs in ~0.7 s where nuSQuIDS needs ~10 minutes to reach a tolerance at which
its output is even a probability.

## Should NuOscProbExact be in this figure?

It has no route to this object: it returns `P(E, L)` with full coherence. A
number can be *derived* from it by averaging over many oscillation cycles, and
asymptotically that converges to the same limit -- but only if two things hold,
and both are measurable rather than assumable.

**Gate, to run before committing to the figure.** Compute Magnus's averaged
probability and, beside it, a mean of NuOscProbExact's coherent output over an
increasing number of samples at the same production point. Then ask:

1. **Do they converge to each other**, or plateau at different values? They agree
   only if every eigenvalue pair is fully decohered over the averaging window. The
   `avgprob` module deliberately keeps cross terms for pairs that stay *coherent*,
   and a wide-window numerical mean destroys exactly those. At solar energies over
   1 AU the pairs are expected to be fully decohered, so the blocks should all be
   singletons -- but "expected" is the word that has cost this project the most
   time, and the module exists because that assumption is not always right.
2. **Does the averaging window matter?** If the converged value depends on how
   wide a window is averaged, the two codes are not computing the same thing and
   the comparison is a category error dressed as a measurement.

**If the gate passes**, include NuOscProbExact with the number of averaging
samples as its dial -- the direct analogue of Magnus's `rtol`. More samples cost
more time and give a better mean, which is a genuine speed-against-accuracy
curve of exactly Figure 11's kind. The caption and the text must then say
plainly that the averaged probability was computed *by us* from its output, and
how.

**If the gate fails, drop it**, and say why in one sentence: the figure then
carries a single code, like the bottom row of Figure 11, and the honest statement
is that no closed-form competitor computes this object at all.

**A caveat that survives either outcome.** Even when the gate passes, the x-axis
for the derived curve measures *our sampling choice*, not NuOscProbExact's
efficiency at a task it was built for. The caption should not imply otherwise.

## The reference

This is the part with no precedent in the existing figures, and it needs a
decision.

A DOP853 or mpmath integration of the coherent evolution is **not** a reference
for this quantity: it returns an oscillating probability, and the averaged one is
its decohered limit, not a converged value of it.

Two candidates:

1. **A high-precision evaluation of the same closed form** -- the eigenbasis
   decomposition at production and detection, and the crossing probabilities,
   carried in mpmath at high precision. This scores the implementation, not the
   physics: it answers "is the closed form evaluated accurately", which is the
   right question for a speed-against-accuracy panel, and it mirrors what
   `prem_chord_reference.json` does for the PREM figure.
2. **A converged numerical average of a high-precision coherent solution** --
   independent of the closed form, so it also tests whether the closed form is
   the right limit. Far more expensive, and it inherits the window question above.

**Recommendation: (1) as the reference, with (2) run once as a cross-check at a
single configuration** and reported in the text rather than plotted. That gives
the panel a reference of the same kind as the other figures while still asking,
once, whether the closed form is the object it claims to be.

## Structure

Five rows, one column: the BS2005-AGS,OP profile at the top, then the averaged
probability at 2, 3, 3+1 and 3+2 flavours. Figure 11's plotting cell is already
parameterised by `(col_key, axes, bench)`, so a one-column variant is a small
edit rather than a new figure. `RTOL_LABEL_OFFSETS` needs one new key.

## Cost

Unknown and probably modest on the Magnus side -- the averaged probability is
closed-form, and the whole point of the figure is that it is cheap. The expensive
parts are the reference (mpmath at high precision, eight cases) and, if the gate
passes, NuOscProbExact's sampled means at the densest settings.

Measure both codes in one session, as Figure 11 now does; this machine has been
measured drifting 12-20% on overhead-bound work.
