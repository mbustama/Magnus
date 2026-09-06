# The 3.49 x 2.1 != 3.66 discrepancy, resolved

**Speedup ratios measured against a shared endpoint do not multiply; time saved adds.
The 2x2 factorial shows the two kernels' savings are additive to within noise, and the
"2.1x" commutator figure is reproduced exactly as the step from compose-on to both-on.**

## The measurement

All four kernel on/off combinations in one process, PREM chord, order 4, fixed
544->1088, marginal us/slab/E, control ratio 1.02 (`scratchpad/survey/factorial.py`):

| d | none  | compose only | comm only | both  |
|---|-------|--------------|-----------|-------|
| 2 | 0.979 | 0.537 (1.82x)| 0.610 (1.61x) | 0.186 (5.26x) |
| 3 | 1.451 | 0.744 (1.95x)| 1.089 (1.33x) | 0.360 (4.03x) |

Additivity: at d=3 the composition kernel saves 0.707, the commutator kernel saves
0.362, and none - both = 1.091 vs 0.707 + 0.362 = 1.069 -- additive to 2%. Same at d=2
(0.443 + 0.369 = 0.812 vs 0.793). There is no overlap between the gains; nothing is
mysterious about the times themselves.

## Why the product predicted wrong

With T = c + k + r (composition cost, commutator cost, everything else):

- the commutator's "about 2.1x" was evidently measured **with the composition kernel
  already on**: (c... = compose-on)/(both) = 0.744/0.360 = **2.07x** -- the match is
  exact. It is (k + r)/r.
- the composition's 3.49x is consistent with being measured **with the commutator
  kernel on** (1.089/0.360 = 3.02x here; his 3.49 vs my 3.02 is regime drift -- his
  runs used different slab counts, refinement on in places, and a 1.5x control drift).
  It is (c + r)/r.
- their product is (c+r)(k+r)/r^2, but the true combined figure is (c+k+r)/r. The
  ratio between them is 1 + ck/(rT): the product **double-counts the shared small
  endpoint r**. At d=3, ck/(rT) = (0.707*0.362)/(0.360*1.451) = 0.49, so the product
  overpredicts by ~1.5x -- which is the whole discrepancy (7.3 predicted vs ~4 actual).
- at d=2 it "roughly matched" by accident: his individual d=2 figures came from regimes
  whose drift happened to cancel against the double-count. The factorial shows d=2 is
  just as additive.

## The rule worth keeping

When several changes each remove an additive share of the same loop, report **time
saved per slab**, not speedup ratios -- savings compose by addition regardless of what
else has landed, and ratios compose only through the shrinking remainder. (This is also
why the combined stack in reports 01-03 is quoted in us/slab/E throughout.)
