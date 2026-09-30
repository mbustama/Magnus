# Measurement records

The scripts and results behind numbers that the package documents (docstrings, `CHANGELOG.md`,
`docs/source/`, `docs/dev/`) or that the paper cites, one folder per piece of work.  They are
kept as they were run.  Each ran against the commit its README names, and most point at the
scratch paths of the session that ran them, so they record the protocol rather than run as
they are; the READMEs say which are worth rerunning and how.

Harnesses maintained as tools live in `tools/` instead: `tools/auto_tight/` holds #120's
measurement (runner, case list, DOP853 references).

| Folder | What it backs |
|---|---|
| `issue070_auto_ladder_route/` | `AUTO_LADDER_MAX_PHASE`, `AUTO_LADDER_MAX_FLOOR_FRACTION`, the #70 entries of `CHANGELOG.md` and `DECISION_DISPATCH_ORDER.md` |
| `notebook_route_audit/` | Which engine every notebook call reaches; used for #70, #122 and #120 |
| `issue089_hybrid_patch_tolerance/` | The hybrid patch-tolerance fix (PR #89) |
| `issue091_093_adoption/` | The re-measured #91 series and #93 tolerance series in `notebooks/external_*.json` (PR #109) |
| `issue092_paper_cache_recompute/` | The recompute of `notebooks/paper_figure_cache.json` (PR #119) |
| `issue120_auto_tight/` | `AUTO_LADDER_TIGHT_MAX_PHASE` and the #120 text in `CHANGELOG.md`, `engines.rst` and `DECISION_DISPATCH_ORDER.md` (PR #128) |
| `issue122_clamp_refusal/` | The clamp table in `MIN_EFFECTIVE_REFINEMENT`'s docstring and the #122 entry of `CHANGELOG.md` (PR #123) |
| `issue160_audit/` | The #160 checklist run case by case; the reproducer table and "Closes" line of the bundled #155/#160 pull request |
| `issue161_rtol_gap/` | The measured true-error-to-tolerance ratios on smooth PREM chords in `diagnostics.rst` ("What `rtol` and `atol` actually control") |
| `issue161_undeclared_features/` | How often undeclared jumps and narrow spikes ended the ladder on a wrong answer, and the 140-call benchmark (values and timing) behind the #161 entry of `CHANGELOG.md` |
| `issue165_pseudo_dirac_resolution/` | The accuracy of `average=True` on a pseudo-Dirac pair, generic route against `osc_prob_pseudo_dirac_vacuum`, in their docstrings, `CHANGELOG.md` and `averaged_probability.rst` |
| `listing1_timing/` | The paper's timing sentence for Listing 1 |

The plan these came from is `docs/dev/HANDOVER_AUTO_TIGHT.md`.  Not kept, and archived outside
the repository with the rest of `resources/handover_auto_tight/`: test and docs-build logs,
process files, drafts of docstrings now in the code, the sessions' transcripts, and a 16 MB
recomputed cache that PR #119 merged.  The folder is excluded from `ruff` in `pyproject.toml`,
like `docs/dev/overhead_survey/prototypes/`, so that the scripts stay as they were run.
