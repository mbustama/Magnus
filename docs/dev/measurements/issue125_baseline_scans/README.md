# Issue #125: baseline scans under `strategy='auto'` at tight tolerances

Source of `HYBRID_YIELDS_TO_CUMULATIVE_MIN_POINTS_TIGHT` and of the #125 text in `engines.rst`,
`performance.rst`, `DECISION_DISPATCH_ORDER.md` and `CHANGELOG.md`.  Measured on `main` at
7e2f13a; `auto` is `main`, `prop` is `main` with the threshold set to 2 (baseline scans) or
`strategy='magnus'` (energy scans).  Every error is against `solve_ivp`/DOP853 at
`rtol=1e-13`, `atol=1e-15`, with a second run at 1e-12/1e-14 as the reference's own spread.

- `b125.py`: the battery.  #120's workload families (exponential profiles 25 to 25 000 km at two
  to five flavors, NSI, LIV, a multi-resonance profile, coherent solar chords), as baseline scans
  of 2, 4 and 7 points at `rtol` = 1e-7, 1e-9, 1e-12 (`atol` a hundredth), orders 4 and 8; and
  energy scans of 16 energies.  `python b125.py run OUT.jsonl NPROC`, `python b125.py score OUT.jsonl`.
  `b125_battery.jsonl` is its output: 96 baseline-scan records and 7 energy-scan records.  The
  coherent solar chords to the full radius (phase 1e5 to 1e6 rad) were stopped, not finished:
  their DOP853 references do not end in useful time, and phase-averaged solar requests
  (`average=True`), the realistic setting, never reach this choice.
- `b125p.py`: the physical population of `docs/dev/adversarial_batteries/physical_profiles.py`
  (tabulated, BS05, supernova shock w = 1e-2 and 1e-6, turbulence C* = 0.1, an Earth crust with
  undeclared steps), lowest energy of each family (45 MeV for the supernova families, whose
  phase at 15 MeV is above 5e4 rad), 2 and 7 baselines, `rtol` = 1e-7 and 1e-12, order 4.
  `b125_physical.jsonl` is its output; score it with `b125.py score`.

Result (`b125.py score`): over 110 baseline scans the cumulative scan took 0.06 to 0.18 of the
hybrid's time at the median, at most 1.08 (where both sides already ran it), and made no silent
miss the hybrid did not also make.  It added 3 to 18 `ToleranceNotAchievedWarning`s per setting
of the 96, all on answers inside the request.  The energy-scan half (`strategy='magnus'` against
`auto`) was measured on 7 non-solar workloads only, too few to touch `AUTO_LADDER_TIGHT_MAX_PHASE`,
which stays.
