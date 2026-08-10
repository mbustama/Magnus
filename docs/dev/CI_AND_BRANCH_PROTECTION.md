# CI and branch protection

**Written:** 2026-08-10, when `main` was first protected.

`main` is protected by a repository **ruleset** — not classic branch protection, which is a
different and older API. GitHub → Settings → Rules, or:

```bash
gh api repos/mbustama/Magnus/rulesets/20617286
```

| rule | effect |
|---|---|
| `deletion` | `main` cannot be deleted |
| `non_fast_forward` | no force-pushes to `main` |
| `pull_request` | changes reach `main` through a PR; **0 approvals required** |
| `required_status_checks` | the six below must pass |

Required: `pytest (3.10)`, `pytest (3.11)`, `pytest (3.12)`, `Run Ruff Linter`,
`CLI reference is up to date`, `Documentation builds`.

**Zero required approvals is deliberate**, not an oversight: a solo maintainer who requires an
approval cannot merge their own work at all. What the rule buys here is the *PR itself* — every
change to `main` gets a commit range, a diff and a set of checks — not a second pair of eyes
there is nobody to provide.

**The admin role bypasses always** (`bypass_actors: RepositoryRole 5`), so you can never lock
yourself out of your own repository. Confirm with `current_user_can_bypass` in the API response
above. Use it sparingly; it is an escape hatch, not a workflow.

Mirrors `mbustama/NuOscProbExact`, ruleset 20558423, which is worth diffing against if either
drifts.

## Why `Coverage` is not required

It takes ~22 minutes and would set the floor on every PR's merge latency, while `pytest` already
runs the same tests on three Python versions. It still runs, and it is still worth reading; it is
just not worth blocking on. NuOscProbExact omits it for the same reason.

`build-and-deploy` (Pages) and `update-pip-graph` must **never** be required: they run only on
`main` or on a schedule, so they never report on a PR head, and a required check that never
reports blocks the PR forever.

## The trap: `paths:` on a trigger and required checks do not mix

**A `paths:` filter on a workflow trigger does not skip the job. It stops the workflow from
existing** — no run, and therefore no check run on the commit at all. Verified on this
repository's own history: commits `7cf2a63` and `f8bfdd3`, which touch only
`docs/dev/HANDOVER_OVERHEAD.md`, produced **zero** notebook check runs, while `bb1e2fd`, which
touched `src/`, produced seven.

That is harmless until one of those contexts is *required*, at which point a pull request
touching only `docs/` or `README.md` waits forever for a check that will never be reported —
**pending, not failing**, so the PR gives the reader nothing to act on and no reason why.

`notebooks.yml` filtered this way until PR #38. It now filters *inside* the workflow: a
`Detect what changed` job answers the same four-path question, the two notebook jobs run on
every commit, and only their expensive **steps** are skipped. The contexts are therefore always
reported and can be required. The cost is a few seconds of runner start-up on commits that
touch no notebook.

**So: do not add a trigger-level `paths:` filter to any workflow whose checks are required, and
do not "simplify" `notebooks.yml` back into one.** If a job must be conditional, make the *job*
conditional or the *steps* conditional — both still report.

The gate fails safe in both directions, and both directions matter:

* Anything the detector cannot resolve — a new branch, a force push, a base commit that is not
  fetched — reports `true`, and everything runs.
* The step guards test `!= 'false'` rather than `== 'true'`, so if the detector fails outright
  and its output is absent, the notebooks still execute.

A gate like this may err towards doing too much. It must never silently do nothing.

## Adding a required check

Add the context only once you have seen it report on a commit that does **not** exercise it —
that is the case that blocks, and the case a normal PR will not reveal. The order is: make the
workflow report unconditionally, watch a PR that skips the work, *then* require it.
