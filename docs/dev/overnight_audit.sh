#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
#
# Unattended pre-publish verification, in two phases.
#
# ---------------------------------------------------------------------------
# WHY TWO PHASES, WHICH IS THE WHOLE POINT OF THIS SCRIPT
# ---------------------------------------------------------------------------
# Phase 1 is correctness: pass/fail work that does not care what else the
# machine is doing.  Phase 2 is cost: benchmarks whose *only* output is a
# timing, which are worthless if anything else is running.
#
# This is not hypothetical.  P4 was once run for 1 h 41 m alongside a
# `pytest -n auto` across twelve cores and the result had to be thrown away --
# and the numbers looked entirely plausible while being meaningless, which is
# the dangerous part.  So phase 2 refuses to start until the machine is quiet,
# and runs strictly one job at a time.
#
# ---------------------------------------------------------------------------
# SAFEGUARDS, because this runs while nobody is watching
# ---------------------------------------------------------------------------
#   * every job is wrapped in run_capped.sh, which sets MemoryMax and
#     MemorySwapMax=0 -- a runaway job is killed rather than left to thrash the
#     machine into a swap death that needs a reboot
#   * every job has a hard timeout, so a hang costs one slot rather than the night
#   * a failing job is recorded and the run continues; one broken thing must not
#     cost the other six
#   * preflight refuses to start on a dirty tree, on low disk, or on low memory
#   * logs and the summary go outside the repository, so nothing here can trip
#     the file-tree guard or end up in a commit
#
# Usage:  docs/dev/overnight_audit.sh [--phase1-only|--phase2-only] [--dry-run]

set -uo pipefail        # NOT -e: a failing job must not abandon the rest

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CAPPED="$REPO/docs/dev/run_capped.sh"
STAMP="$(date +%Y%m%d-%H%M%S)"
LOGDIR="$HOME/magnus_overnight/$STAMP"
SUMMARY="$LOGDIR/SUMMARY.txt"

PHASE1=1; PHASE2=1; DRYRUN=0
for arg in "$@"; do
    case "$arg" in
        --phase1-only) PHASE2=0 ;;
        --phase2-only) PHASE1=0 ;;
        --dry-run)     DRYRUN=1 ;;
        *) echo "usage: $0 [--phase1-only|--phase2-only] [--dry-run]" >&2; exit 2 ;;
    esac
done

# Caps chosen against a 15 GB machine with ~7 GB free.  The suite under
# `-n auto` is the hungriest; the batteries are single-process and modest.
CAP_HEAVY=5G            # notebook build, coverage, full suite
CAP_LIGHT=3G            # batteries, calibration, clean-room install

# Phase 2 will not start above this 1-minute load average.  1 means "one core
# busy on a twelve-core box", which is already enough to move a benchmark;
# raise it only if you accept noisier timings.  Override: QUIET_LOAD=2 ...
QUIET_LOAD="${QUIET_LOAD:-1}"

mkdir -p "$LOGDIR"

say() { printf '%s\n' "$*" | tee -a "$SUMMARY"; }

# ---------------------------------------------------------------------------
# preflight
# ---------------------------------------------------------------------------
preflight() {
    local fail=0
    say "=== preflight ==="
    say "repository : $REPO"
    say "commit     : $(cd "$REPO" && git rev-parse --short HEAD 2>/dev/null || echo '?')"
    say "branch     : $(cd "$REPO" && git rev-parse --abbrev-ref HEAD 2>/dev/null || echo '?')"

    if [ -n "$(cd "$REPO" && git status --porcelain 2>/dev/null)" ]; then
        say "REFUSING: the working tree is dirty.  Results have to correspond to a"
        say "          commit, or a failure cannot be attributed to anything."
        fail=1
    fi

    local disk_gb
    disk_gb=$(df -BG --output=avail "$REPO" | tail -1 | tr -dc '0-9')
    say "disk free  : ${disk_gb} GB"
    # the notebook build rewrites every .ipynb with outputs, and a render writes
    # ~224 MB of raw GIFs before they are shrunk
    if [ "${disk_gb:-0}" -lt 5 ]; then
        say "REFUSING: under 5 GB free.  The notebook build and the renders need room,"
        say "          and running out mid-way leaves half-written notebooks."
        fail=1
    fi

    local mem_gb
    mem_gb=$(free -g | awk '/^Mem:/ {print $7}')
    say "mem avail  : ${mem_gb} GB"
    if [ "${mem_gb:-0}" -lt 4 ]; then
        say "REFUSING: under 4 GB available.  Every cap below assumes more headroom"
        say "          than that, and MemorySwapMax=0 turns a squeeze into a kill."
        fail=1
    fi

    if [ ! -x "$CAPPED" ]; then
        say "REFUSING: $CAPPED is missing or not executable."
        fail=1
    fi
    say ""
    return $fail
}

# ---------------------------------------------------------------------------
# one job: capped, timed out, logged, non-fatal on failure
# ---------------------------------------------------------------------------
run_job() {
    local name="$1" cap="$2" timeout_s="$3"; shift 3
    local log="$LOGDIR/${name}.log"
    local started elapsed rc

    say "--- $name  (cap $cap, timeout $((timeout_s/60)) min)"
    if [ "$DRYRUN" -eq 1 ]; then
        say "    DRY RUN: $*"
        return 0
    fi

    started=$(date +%s)
    ( cd "$REPO" && timeout --signal=TERM --kill-after=60 "$timeout_s" \
        "$CAPPED" -m "$cap" -- "$@" ) >"$log" 2>&1
    rc=$?
    elapsed=$(( $(date +%s) - started ))

    case $rc in
        0)   say "    OK        ($((elapsed/60)) min)  -> ${name}.log" ;;
        124) say "    TIMEOUT   ($((elapsed/60)) min)  -> ${name}.log" ;;
        137) say "    KILLED    ($((elapsed/60)) min)  -- hit the memory cap; see ${name}.log" ;;
        *)   say "    FAILED rc=$rc ($((elapsed/60)) min)  -> ${name}.log"
             tail -5 "$log" | sed 's/^/      | /' | tee -a "$SUMMARY" >/dev/null ;;
    esac
    return 0
}

# ---------------------------------------------------------------------------
# phase 2 gate: the machine has to be quiet, or the timings are fiction
# ---------------------------------------------------------------------------
wait_for_quiet() {
    local tries=0 load
    if [ "$DRYRUN" -eq 1 ]; then
        say "DRY RUN: would wait here for the machine to go quiet"
        return 0
    fi
    while [ $tries -lt 60 ]; do
        load=$(awk '{print int($1+0.5)}' /proc/loadavg)
        if [ "$load" -le "$QUIET_LOAD" ] && ! pgrep -f 'pytest|make_notebooks|sphinx' >/dev/null; then
            say "machine quiet (load ${load}); starting the cost measurements"
            return 0
        fi
        [ $tries -eq 0 ] && say "waiting for the machine to go quiet (load ${load})..."
        sleep 60; tries=$((tries+1))
    done
    say "STILL BUSY after an hour -- skipping phase 2 rather than recording"
    say "timings that cannot be trusted.  Re-run with --phase2-only when idle."
    return 1
}

# ===========================================================================
say "Magnus overnight verification -- $STAMP"
say "logs: $LOGDIR"
say ""

preflight || { say "preflight failed; nothing run."; exit 1; }

if [ "$PHASE1" -eq 1 ]; then
    say "=== PHASE 1: correctness (may share the machine) ==="

    # Fast things first, so a broken one is visible early rather than at 4 a.m.
    run_job 01-sev-tol-calibration "$CAP_LIGHT" 1800 \
        python3 docs/dev/calibrate_sev_tol.py --bases 8

    run_job 02-clean-room-install "$CAP_LIGHT" 2400 \
        bash -c 'set -e
            rm -rf /tmp/magnus_cleanroom && mkdir -p /tmp/magnus_cleanroom
            python3 -m build --outdir /tmp/magnus_cleanroom/dist
            python3 -m venv /tmp/magnus_cleanroom/venv
            /tmp/magnus_cleanroom/venv/bin/pip -q install --upgrade pip
            /tmp/magnus_cleanroom/venv/bin/pip -q install "/tmp/magnus_cleanroom/dist"/*.whl pytest pytest-xdist
            # the suite against the INSTALLED package, not the source tree: this is
            # what catches a data file that never made it into the wheel
            cd /tmp && /tmp/magnus_cleanroom/venv/bin/python -m pytest "'"$REPO"'/tests" -q -n auto'

    run_job 03-coverage "$CAP_HEAVY" 3600 \
        python3 -m pytest tests/ -q -n auto --cov=magnus --cov-report=term

    # The long pole.  Executes all 27 notebooks, which is what CI does and what
    # the two results-changing defaults (NuFIT 6.1, Earth Y_e) most need.
    run_job 04-notebooks-all "$CAP_HEAVY" 21600 \
        python3 notebooks/make_notebooks.py

    run_job 05-docs-clean "$CAP_LIGHT" 3600 \
        bash -c 'cd docs && rm -rf build source/api &&
                 python3 -m sphinx -b html -n -W --keep-going source build'

    say ""
fi

if [ "$PHASE2" -eq 1 ]; then
    say "=== PHASE 2: cost (strictly alone) ==="
    if wait_for_quiet; then
        # P4: the seam-cost battery, recorded as never having been run to
        # completion.  ~3 h, and it carries its own interleaved control.
        run_job 10-P4-seam-cost "$CAP_LIGHT" 18000 \
            bash -c 'cd docs/dev/adversarial_batteries && python3 physical_battery.py seam_cost'

        # The PREM-touching batteries.  Their recorded findings were all
        # established with a uniform Y_e = 0.5, which no longer holds.
        run_job 11-physical-profiles "$CAP_LIGHT" 7200 \
            bash -c 'cd docs/dev/adversarial_batteries && python3 physical_profiles.py'

        run_job 12-validate-physical "$CAP_LIGHT" 7200 \
            bash -c 'cd docs/dev/adversarial_batteries && python3 validate_physical.py'

        # bitident.py dumps its raw probabilities to a path it takes as argv[1];
        # called bare it dies on IndexError before running anything.  The .npz goes
        # beside the log, i.e. outside the repository, for the same reason the logs
        # do: nothing this script writes may reach the file-tree guard or a commit.
        run_job 13-bitident "$CAP_LIGHT" 3600 \
            bash -c "cd docs/dev/adversarial_batteries && python3 bitident.py '$LOGDIR/13-bitident.npz'"
    fi
    say ""
fi

say "=== done: $(date +%H:%M:%S) ==="
say ""
say "Read $SUMMARY first; each job's full output is beside it."
say ""
say "Two things this script cannot do for you:"
say "  * it re-executes the notebooks but does not READ them.  A notebook can"
say "    run clean and still assert something in prose that its new figures"
say "    contradict -- the Earth and NuFIT changes both move numbers."
say "  * it does not re-render anim_earth.gif, which still shows the old"
say "    composition.  That needs RENDER = True in notebook 27."
