# Issue #70: `strategy='auto'` hands a moderate phase to the ladder (PR #72)

Run during issue #70, before PR #72 merged it.  Paths inside point into that session's scratch.

- **The campaign:** `campaign.py` writes `campaign.jsonl`, 87 records over 20 smooth workloads
  (exponential profiles from 25 to 250 000 km, 2 to 5 flavors, NSI, LIV) with DOP853 references at
  1e-12: the hybrid strategy against the ladder at a tenth of the tolerance.  `score.py` scores
  it; `crossover.py` prices both routes across the phase, point by point.  These back
  `AUTO_LADDER_MAX_PHASE`'s docstring and the #70 entries of `CHANGELOG.md` and
  `DECISION_DISPATCH_ORDER.md`.
- **Route validation:** `validate.py` and `validate_sun.py` write `validate.jsonl`;
  `validate_v1_floor_oldphase.jsonl` is the first version, with the old phase measure.
- **Phase measures:** `measures.py` writes `measures.jsonl`: the norm of the integrated
  Hamiltonian against the spread of its eigenvalues, which the dispatcher now uses.
- **Which engine answers under the rule:** `routes.py`.
- **Listing 1:** `check_l1.py`, `rungs.py` and `l1_ref.npz`; `engines_both.py` (every engine at
  both tolerances); `order_tight.py` (orders 4, 6 and 8 at 1e-12); `worst_dop.py`;
  `ladder_check.py`, which writes `ladder_check_*.json` (the order-8 listing against DOP853 at 26
  energies per curve); `fig1_vs_listing.py` (the listing against the cached Fig. 1 data).
- **The Sun and the patch fix:** `sun_floor.py` (two flavors, 10 MeV over 0.9 R_sun, behind
  `AUTO_LADDER_MAX_FLOOR_FRACTION`), `sun10.py`, `sun_auto.py`, `sun_patch.py`; `patch_probe.py`,
  the root-cause probe of PR #89, varying only `patch_atol`.
- **The notebook route audit's record:** `route_audit.log` (four notebooks rerouted: 13, 20, 22
  and 27).  The audit is in `../notebook_route_audit/`.
- **History only:** `nb24_shock.py`; `patch_floor.py` and `patch_narrow.py`, source edits since
  applied.

Worth rerunning, with its paths fixed, when the route's thresholds are revisited: `campaign.py`.
