# Development history removed from the user docs

The final docs audit (issue #193) took development history out of the user pages. Most of it
was already in `CHANGELOG.md`, elsewhere in `docs/dev`, or in a constant's docstring. The
passages below were recorded nowhere else, so they are kept here as they were written.

## conventions.rst

> Magνs has been bitten by exactly that: a reversed slab ordering, a doubled antineutrino
> potential sign and a flipped two-flavor mass ordering were all fixed on the same day, and
> each had been silently self-consistent.

## adiabatic_strategy.rst

The complex-step derivative:

> This was an early implementation bug, caught by comparing against a real-valued toy
> Hamiltonian (where the same complex-step formula happened to work by coincidence, since the
> function actually was real-valued there) -- a cautionary example of why this module never
> uses it, on any Hamiltonian, real or complex.

The resonance threshold:

> this was checked directly during development: a fixed `threshold=0.1` gave a 2% error on one
> engineered case, needing `0.01` for better than 1e-4

## methodology.rst

The seed prototype:

> an order-aware target was prototyped and A/B tested over 45 configurations (five cases ×
> three orders × three tolerances). It gave no speed-up, and cost up to 20% on the energy
> scan: the final slab count is set by the refinement loop, not the seed, so starting coarser
> only adds an iteration.

## engines.rst

The tight-tolerance ladder route (issue #120; raw data in `measurements/issue120_auto_tight/`):

> on 300 energies from 3 to 100 MeV over 200 km of an exponential profile (418 rad) at 1e-8,
> 0.1 s where the adiabatic engine took 11 s

## architecture.rst

> When every wrapper declared its own copy of these keywords, the copies drifted: one wrapper
> had a different default tolerance, one lacked `nubar`, one had a different validation bound.

> four functions in this package did exactly that before a check was written for it.

## plotting.rst

> this project has already paid for that: `oscprob`'s keyword chain used to forward unknown
> names down several layers before failing somewhere unrecognizable.

From "Why it exists":

> a nine-keyword `legend` invocation, four `MultipleLocator` assignments … copied from figure
> to figure and varied slightly each time
