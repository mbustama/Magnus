#!/usr/bin/env bash
# Runs a long job under a hard memory cap, so that a runaway dies alone.
#
# WHY THIS EXISTS.  Three separate jobs in one session outgrew this machine, and the third
# took the whole editing session down with it:
#
#   * a notebook-25 benchmark that asked for a single 6.0 GB array and was killed by the
#     kernel 25 minutes in, leaving a zero-byte output file;
#   * a DOP853 reference at an eV-scale sterile splitting, which needed ~5.9e6 oscillations
#     resolved and would have run for about a day;
#   * `make clean html`, which is not just Sphinx -- docs/source carries 22 `jupyter-execute`
#     directives, each spawning a Jupyter kernel that imports Magnus and computes.  Sphinx,
#     AutoAPI and those kernels together exhausted memory, and the kernel's OOM killer chose
#     its victim SYSTEM-WIDE rather than inside the build.
#
# The fix is containment, not politeness: a systemd scope with MemoryMax means the cgroup's
# own OOM killer fires first and kills only what is inside the scope.  Verified rather than
# assumed -- allocating 900 MB under a 512 MB cap kills the child and nothing else.
#
# Single-threaded BLAS is not a micro-optimisation here.  On a 12-core box every kernel
# otherwise spawns a full thread pool, multiplying both the memory and the contention, and
# for these jobs the arithmetic is already batched at the numpy level.
#
# USAGE
#     docs/dev/run_capped.sh [-m 6G] -- <command> [args...]
#
#     docs/dev/run_capped.sh -m 6G -- make -C docs clean html SPHINXOPTS="-n -W --keep-going"
#     docs/dev/run_capped.sh -m 4G -- python notebooks/gen_shock_benchmarks.py > out.json
#
# Pick the cap from what is FREE (`free -g`), not from what is installed, and leave a couple
# of gigabytes for everything else.

set -euo pipefail

CAP=6G
while [ $# -gt 0 ]; do
    case "$1" in
        -m) CAP="$2"; shift 2 ;;
        --) shift; break ;;
        *)  echo "usage: $0 [-m CAP] -- <command>" >&2; exit 2 ;;
    esac
done

if [ $# -eq 0 ]; then
    echo "usage: $0 [-m CAP] -- <command>" >&2
    exit 2
fi

if ! command -v systemd-run >/dev/null 2>&1; then
    # Better to run uncapped than to refuse: say so loudly rather than failing a build.
    echo "run_capped: systemd-run not available; running WITHOUT a memory cap" >&2
    exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        NUMEXPR_NUM_THREADS=1 "$@"
fi

echo "run_capped: MemoryMax=$CAP, swap off, single-threaded BLAS, nice 10" >&2
exec systemd-run --user --scope \
    -p MemoryMax="$CAP" -p MemorySwapMax=0 --nice=10 --quiet -- \
    env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        NUMEXPR_NUM_THREADS=1 "$@"
