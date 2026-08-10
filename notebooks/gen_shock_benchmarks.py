# -*- coding: utf-8 -*-
r"""gen_shock_benchmarks.py

Speed against accuracy for Mag(nu)s and NuOscProbExact on the supernova shock ray of
notebook 14, plus a probability-against-energy panel.  Writes
``external_shock_benchmarks.json``, which notebook 25 reads; neither external code is
needed to *run* that notebook.

    python notebooks/gen_shock_benchmarks.py > notebooks/external_shock_benchmarks.json

THE PHYSICS IS NOT COPIED.  Notebook 14 already defines this profile, the ray, the
energy and the sampled baselines, and it has a frozen DOP853 reference to match.
Transcribing any of that would create a second definition that can drift from the first,
so this script *executes notebook 14's own cells* and reads the definitions out of the
resulting namespace -- the same device ``make_shock_reference.py`` uses, and for the same
reason.  The referee is that notebook's frozen ``solve_ivp``/DOP853 solution, loaded
through its own fingerprint guard, so a profile change invalidates it loudly.

TWO SETUP TRAPS WERE HIT WHILE WRITING THIS, both of the kind that converge cleanly to
the wrong answer.

  Mag(nu)s must be driven CUMULATIVELY.  The case asks for the probability at 61 points
  along one ray.  Driven point by point, Mag(nu)s re-propagates the whole ray for each
  one: 3.27 s and an error of 4.7e-04 that does not move with ``rtol``, because the
  refinement ladder is being asked for a tolerance it cannot reach that way.  With
  ``cumulative=True`` the same request costs 0.24 s and lands at 8.1e-06.  A comparison
  against the point-by-point route would have been a comparison against the wrong driver.

  NuOscProbExact must be given enough slabs INSIDE THE FRONT.  The two fronts are 0.07 km
  of a 70 000 km ray, so slabs allocated in proportion to length leave each front with a
  single sample and the code appears to floor at 5e-07 no matter how fine the rest gets.
  That floor is this script's allocation, not NuOscProbExact's limit: raising the
  per-front minimum to 1/1000 of the budget takes the same code to 1.4e-09.  The
  allocation below is therefore a choice made on its behalf, and it is stated in the
  notebook rather than hidden here.

HOW EACH CODE IS DRIVEN, and why that is the way its authors intend.

  Mag(nu)s: ``cumulative=True`` with the fronts declared through ``t_breakpoints``, one
  call for all 61 baselines.  Dial is ``n_slabs``.

  NuOscProbExact: ``evolution_operator_3nu_slabs`` accumulated.  Its batched route returns
  the operator at the END of a slab chain and this case wants it at 61 points WITHIN one,
  so the ray is cut at the declared fronts and at every target, each leg solved, and the
  operators composed.  The 60 intervals between consecutive targets are equal by
  construction (``Ls`` is a linspace), so subdivided equally they share one set of widths
  and go in ONE batched call -- looping them would time the loop rather than the code.
  Dial is the total slab count.
"""

__author__ = "Mauricio Bustamante"
__email__ = "mbustamante@gmail.com"


import json
import pathlib
import platform
import sys
import time
import warnings

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'src'))
sys.path.insert(0, str(HERE))

import magnus.globaldefs as gd                                    # noqa: E402
import magnus.oscprob as oscprob                                  # noqa: E402

NOTEBOOK = '14_magnus_supernova_shock.ipynb'
WIDTHS = (1e-6, 1e-3)
MG_DIALS = (2000, 8000, 32000, 128000, 512000)
NPE_DIALS = (2048, 8192, 32768, 131072, 524288)
N_ENERGY = 120
E_LO_MEV, E_HI_MEV = 5.0, 60.0
MG_ENERGY_SLABS = 128000
NPE_ENERGY_SLABS = 524288
# Energies at which the code-against-code residual is itself refereed by DOP853.  Three,
# because each costs about a minute: the point is to say WHOSE the residual is, not to
# draw a third curve.
E_REFEREED_MEV = (8.0, 15.0, 40.0)

# THE MEMORY GUARD, and it is not decoration -- an earlier run of this script was killed
# by the kernel partway through the energy panel, silently, after twenty-five minutes and
# with an empty output file.  A batch of Hamiltonians is a complex128 array of shape
# (n_energies, n_slabs, 3, 3), which is 144 bytes per energy per slab, so the panel's
# natural phrasing asks for
#
#     Magnus            120 energies x 128 000 slabs  =  2.1 GB
#     NuOscProbExact    120 energies x 372 000 slabs  =  6.0 GB   (its longest segment)
#
# as SINGLE allocations, before either engine's own temporaries, on a machine with about
# 8 GB free.  Both are chunked below instead.  Chunking is exact rather than approximate:
# the energy axis is independent point by point, and consecutive slab chunks compose by
# multiplying their evolution operators, which is the same product in the same order.
MAX_STACK_BYTES = 192*1024**2
BYTES_PER_SLAB_PER_ENERGY = 9*16                  # 3x3 complex128


def chunk_energies(n_energies, n_slabs):
    """How many energies may share one Hamiltonian stack of `n_slabs` slabs."""
    per_energy = max(1, n_slabs*BYTES_PER_SLAB_PER_ENERGY)
    return max(1, min(int(n_energies), int(MAX_STACK_BYTES//per_energy) or 1))


def chunk_slabs(n_energies, n_slabs):
    """How many slabs may share one Hamiltonian stack across `n_energies` energies."""
    per_slab = max(1, int(n_energies)*BYTES_PER_SLAB_PER_ENERGY)
    return max(1, min(int(n_slabs), int(MAX_STACK_BYTES//per_slab) or 1))


def have(name):
    try:
        __import__(name)
        return True
    except Exception:
        return False


HAVE_NPE = have('slabs')


def notebook_namespace():
    r"""Executes notebook 14's definition cells and returns their namespace.

    Stops at ``measure``, which is the first cell that wants a result rather than a
    definition -- exactly where ``make_shock_reference.py`` stops, and for the same
    reason: everything before it is the configuration, and there is one copy of it.
    """
    import contextlib

    import matplotlib
    matplotlib.use('Agg')                         # the cells import pyplot
    import make_notebooks

    ns = {'__name__': '__notebook__'}
    # The notebook's cells print -- they are written to be read, not imported -- and this
    # script's stdout IS the JSON file.  One stray line ("ray = 4729 oscillation
    # lengths...") made the output unparseable while every number in it was correct, so
    # the cells' stdout goes to stderr with the rest of the progress reporting.
    with contextlib.redirect_stdout(sys.stderr):
        for cell in make_notebooks.books[NOTEBOOK].cells:
            if cell.cell_type != 'code':
                continue
            if 'def measure(' in cell.source:
                break
            exec(compile(cell.source, '<%s>' % NOTEBOOK, 'exec'), ns)
    missing = [n for n in ('sn_shock_ne', 'make_H', 'shock_breakpoints',
                           'frozen_reference', 'L0', 'L1', 'Ls', 'ENERGY', 'params3')
               if n not in ns]
    if missing:
        raise SystemExit('%s no longer defines %s; it has been restructured and this '
                         'script needs updating' % (NOTEBOOK, ', '.join(missing)))
    return ns


NS = notebook_namespace()
L0, L1 = NS['L0'], NS['L1']
Ls = np.asarray(NS['Ls'], dtype=float)
ENERGY, PARAMS3 = NS['ENERGY'], NS['params3']


# ------------------------------------------------------------------------- timing
def timed(call, repeat=3, min_block=0.05, budget=6.0):
    """Best of `repeat` autoranged blocks, first pass discarded.

    The discard is not decoration: the first Mag(nu)s call of a session pays ~0.7 s to
    compile the numba kernel, which a user pays once and not once per call.  `budget`
    caps the repeats for dials where a single call already runs for seconds.
    """
    t0 = time.perf_counter()
    call()
    first = time.perf_counter() - t0
    if first > budget:
        return first
    reps = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        el = time.perf_counter() - t0
        if el >= min_block:
            break
        reps *= 2
    best = el/reps
    for _ in range(repeat - 1):
        if best*reps > budget:
            break
        t0 = time.perf_counter()
        for _ in range(reps):
            call()
        best = min(best, (time.perf_counter() - t0)/reps)
    return best


def control():
    a = np.random.default_rng(0).normal(size=(180, 180))
    return a @ a


# ------------------------------------------------------------------------- drivers
def magnus_along_ray(width_frac, n_slabs):
    """P(L0 -> L) at every L in Ls, in one cumulative call."""
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return np.asarray(oscprob.osc_prob_matter_std_potential(
            3, NS['sn_shock_ne'](width_frac), ENERGY, Ls, PARAMS3, L0=L0,
            density_is_of_number_of_electrons=True,
            t_breakpoints=NS['shock_breakpoints'](width_frac),
            cumulative=True, n_slabs=n_slabs, max_n_slabs=4*n_slabs,
            rtol=1.0e-12, atol=1.0e-14)).reshape(len(Ls), 3, 3)


def npe_along_ray(width_frac, n_slabs):
    """The same 61 probabilities, by composing NuOscProbExact slab-chain operators."""
    import slabs as npe

    H = NS['make_H'](NS['sn_shock_ne'](width_frac))
    bps = np.asarray(NS['shock_breakpoints'](width_frac), dtype=float)
    h = (L1 - L0)/float(n_slabs)
    # Never fewer than a thousandth of the budget in any declared segment: the fronts
    # are 0.07 km of 70 000, and proportional allocation alone leaves them unresolved.
    floor = max(8, n_slabs//1000)

    edges = np.unique(np.concatenate([bps[bps < Ls[0]], [L0, Ls[0]]]))
    edges = edges[(edges >= L0) & (edges <= Ls[0])]

    used = 0
    U = np.eye(3, dtype=complex)
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = max(floor, int(round((hi - lo)/h)))
        e = np.linspace(lo, hi, m + 1)
        U = np.asarray(npe.evolution_operator_3nu_slabs(
            np.asarray(H(0.5*(e[:-1] + e[1:])), dtype=complex), np.diff(e))) @ U
        used += m

    # The 60 equal intervals between targets: one batched call, then compose.
    d_l = Ls[1] - Ls[0]
    m = max(1, int(round(d_l/h)))
    e = np.linspace(0.0, d_l, m + 1)
    mid = Ls[:-1, None] + 0.5*(e[:-1] + e[1:])[None, :]
    h_batch = np.asarray(H(mid.ravel()), dtype=complex).reshape(len(Ls) - 1, m, 3, 3)
    u_batch = np.asarray(npe.evolution_operator_3nu_slabs(h_batch, np.diff(e)))
    used += (len(Ls) - 1)*m

    ops = [U]
    for j in range(len(Ls) - 1):
        U = u_batch[j] @ U
        ops.append(U)
    ops = np.array(ops)
    return np.swapaxes(ops.real**2 + ops.imag**2, -1, -2), used


# -------------------------------------------------------------- probability vs energy
def h_at_energy(width_frac):
    r"""Notebook 14's own Hamiltonian, re-pointed to an arbitrary energy.

    ``make_H`` bakes in ``ENERGY``, and the energy enters only through the vacuum term,
    so ``H_E = H_ENERGY + h_vac (1/E - 1/ENERGY)`` is exact.  Adding a matrix is used
    rather than recovering V_CC by subtracting two nearly equal scalars, and rather than
    rebuilding the potential here -- which would put a second definition of the physics
    in a second file, which is the defect this script exists to avoid.
    """
    h_fixed = NS['make_H'](NS['sn_shock_ne'](width_frac))
    hvac = np.asarray(NS['hvac3'])

    def h_of(energy):
        shift = hvac*(1.0/energy - 1.0/ENERGY)

        def f(l):
            return np.asarray(h_fixed(l)) + shift
        return f
    return h_of


def dop853_at_energy(width_frac, energy):
    """An adaptive DOP853 pass over the whole ray -- neither code's method."""
    from scipy.integrate import solve_ivp

    h_of = h_at_energy(width_frac)(energy)

    def rhs(l, y):
        return (-1j*np.asarray(h_of(l)) @ y.reshape(3, 3)).ravel()

    sol = solve_ivp(rhs, (L0, L1), np.eye(3, dtype=complex).ravel(),
                    rtol=1.0e-12, atol=1.0e-14, method='DOP853')
    if not sol.success:
        raise SystemExit('DOP853 failed at %.1f MeV: %s'
                         % (energy/gd.UNIT_MEV, sol.message))
    u = sol.y[:, -1].reshape(3, 3)
    return float(abs(u[gd.NUE, gd.NUE])**2)


def npe_vs_energy(width_frac, energies, n_slabs):
    """P_ee at the end of the ray, batched over energies within a memory budget.

    Two nested chunkings, both exact.  Energies are independent, so a chunk of them is
    the same calculation on fewer rows.  Consecutive slab chunks compose by multiplying
    their evolution operators in order, which is the product the unchunked call forms
    anyway -- `evolution_operator_3nu_slabs` is itself a chain of exactly this.
    """
    import slabs as npe

    energies = np.atleast_1d(np.asarray(energies, dtype=float))
    h_of = h_at_energy(width_frac)
    bp = np.asarray(NS['shock_breakpoints'](width_frac), dtype=float)
    h = (L1 - L0)/float(n_slabs)
    floor = max(8, n_slabs//1000)
    edges = np.unique(np.concatenate([bp, [L0, L1]]))
    edges = edges[(edges >= L0) & (edges <= L1)]

    out = np.empty(len(energies))
    e_step = chunk_energies(len(energies), min(n_slabs, 1 << 16))
    for i0 in range(0, len(energies), e_step):
        e_chunk = energies[i0:i0 + e_step]
        u_all = np.broadcast_to(np.eye(3, dtype=complex),
                                (len(e_chunk), 3, 3)).copy()
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = max(floor, int(round((hi - lo)/h)))
            e = np.linspace(lo, hi, m + 1)
            m_step = chunk_slabs(len(e_chunk), m)
            for j0 in range(0, m, m_step):
                sub = e[j0:j0 + m_step + 1]
                mid = 0.5*(sub[:-1] + sub[1:])
                stack = np.stack([np.asarray(h_of(en)(mid), dtype=complex)
                                  for en in e_chunk])           # (n_chunk, m_sub, 3, 3)
                u_all = np.asarray(npe.evolution_operator_3nu_slabs(
                    stack, np.diff(sub))) @ u_all
                del stack
        out[i0:i0 + e_step] = np.abs(u_all[:, gd.NUE, gd.NUE])**2
    return out


def prob_vs_energy(width_frac):
    """P_ee at the end of the ray across energy, from both codes.

    There is no frozen referee off 15 MeV, so the 120-point curve is code against code.
    Agreement is worth showing and is NOT the same claim as accuracy, so three of those
    energies are additionally refereed by DOP853 -- which is what says whose the residual
    is.  Without that the reader would be free to assume it is shared.
    """
    energies = np.linspace(E_LO_MEV, E_HI_MEV, N_ENERGY)*gd.UNIT_MEV
    ne = NS['sn_shock_ne'](width_frac)
    bps = NS['shock_breakpoints'](width_frac)

    # An explicit slab count, because rtol is not this route's dial: on a SINGLE
    # baseline over this ray the refinement ladder stalls near 1e-2 and says so
    # (MagnusConvergenceWarning, ToleranceNotAchievedWarning).  Given the slabs it
    # reaches 1e-8 or better, refereed below.
    e_step = chunk_energies(len(energies), MG_ENERGY_SLABS)
    chunks = []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        t0 = time.perf_counter()
        for i0 in range(0, len(energies), e_step):
            chunks.append(np.asarray(oscprob.osc_prob_matter_std_potential(
                3, ne, energies[i0:i0 + e_step], L1, PARAMS3, L0=L0,
                nu_i=gd.NUE, nu_f=gd.NUE,
                density_is_of_number_of_electrons=True, t_breakpoints=bps,
                n_slabs=MG_ENERGY_SLABS, max_n_slabs=4*MG_ENERGY_SLABS,
                rtol=1.0e-12, atol=1.0e-14)).ravel())
            print('    magnus energies %d/%d' % (min(i0 + e_step, len(energies)),
                                                 len(energies)), file=sys.stderr)
        p_mg = np.concatenate(chunks)
        t_mg = time.perf_counter() - t0

    p_npe, t_npe = None, None
    if HAVE_NPE:
        t0 = time.perf_counter()
        p_npe = npe_vs_energy(width_frac, energies, NPE_ENERGY_SLABS)
        t_npe = time.perf_counter() - t0

    refereed = []
    for e_mev in E_REFEREED_MEV:
        energy = e_mev*gd.UNIT_MEV
        truth = dop853_at_energy(width_frac, energy)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            mg = float(np.asarray(oscprob.osc_prob_matter_std_potential(
                3, ne, energy, L1, PARAMS3, L0=L0, nu_i=gd.NUE, nu_f=gd.NUE,
                density_is_of_number_of_electrons=True, t_breakpoints=bps,
                n_slabs=MG_ENERGY_SLABS, max_n_slabs=4*MG_ENERGY_SLABS,
                rtol=1.0e-12, atol=1.0e-14)))
        row = dict(energy_mev=e_mev, dop853=truth, magnus_error=mg - truth)
        if HAVE_NPE:
            row['npe_error'] = float(
                npe_vs_energy(width_frac, np.array([energy]), NPE_ENERGY_SLABS)[0]
                - truth)
        refereed.append(row)
        print('    referee %.0f MeV: magnus %+.2e  npe %+.2e'
              % (e_mev, row['magnus_error'], row.get('npe_error', float('nan'))),
              file=sys.stderr)

    return dict(width=width_frac,
                energy_ev=[float(x) for x in energies],
                magnus=[float(x) for x in p_mg],
                magnus_seconds=t_mg,
                magnus_n_slabs=MG_ENERGY_SLABS,
                npe=None if p_npe is None else [float(x) for x in p_npe],
                npe_seconds=t_npe,
                npe_n_slabs=NPE_ENERGY_SLABS,
                refereed=refereed)


CHECKPOINT = HERE/'external_shock_benchmarks.partial.json'


def checkpoint(out, phase):
    """Writes what is finished so far, so a kill does not cost the whole run.

    The first attempt at this script was killed by the kernel twenty-five minutes in,
    during the energy panel, and left a zero-byte output file: every measurement before
    the crash was lost even though none of it was at fault.  Phases are recorded here as
    they complete and `--resume` reads them back.
    """
    out['completed_phases'] = sorted(set(out.get('completed_phases', [])) | {phase})
    tmp = CHECKPOINT.with_suffix('.tmp')
    tmp.write_text(json.dumps(out, indent=1))
    tmp.replace(CHECKPOINT)                      # atomic: never a half-written file
    print('  [checkpoint: %s]' % phase, file=sys.stderr)


def main():
    resume = '--resume' in sys.argv
    if resume and CHECKPOINT.exists():
        out = json.loads(CHECKPOINT.read_text())
        print('resuming; phases already done: %s'
              % ', '.join(out.get('completed_phases', [])) or '(none)', file=sys.stderr)
        run_phases(out)
        json.dump(out, sys.stdout, indent=1)
        return

    out = {'note': ('shock ray of notebook 14, refereed by its frozen DOP853 solution; '
                    'produced by notebooks/gen_shock_benchmarks.py'),
           'machine': platform.platform(),
           'python': platform.python_version(),
           'numpy': np.__version__,
           'energy_mev': ENERGY/gd.UNIT_MEV,
           'L0_km': L0/gd.UNIT_KM, 'L1_km': L1/gd.UNIT_KM,
           'n_targets': len(Ls),
           'targets_km': [float(x)/gd.UNIT_KM for x in Ls],
           'control_ratio': None,
           'control_note': ('two interleaved copies of one workload neither code '
                            'touches; 1.00 is the evidence the ratios are readable'),
           'cases': [], 'prob_vs_energy': [], 'completed_phases': []}

    # Interleaved, not two sequential calls: sequential reads machine drift BETWEEN them
    # rather than contention DURING them.
    best = {'a': np.inf, 'b': np.inf}
    for _ in range(9):
        for k in ('a', 'b'):
            best[k] = min(best[k], timed(control, repeat=1))
    out['control_ratio'] = best['a']/best['b']
    print('control ratio %.3f' % out['control_ratio'], file=sys.stderr)
    checkpoint(out, 'control')

    run_phases(out)
    json.dump(out, sys.stdout, indent=1)


def run_phases(out):
    done = set(out.get('completed_phases', []))

    for width in WIDTHS:
        if 'case:%.0e' % width in done:
            continue
        ref = NS['frozen_reference'](width)          # fingerprint-guarded loader
        ref_unitarity = float(np.max(np.abs(ref.sum(axis=2) - 1.0)))

        mg = []
        for n in MG_DIALS:
            P = magnus_along_ray(width, n)
            t = timed(lambda n=n: magnus_along_ray(width, n))
            mg.append(dict(label=str(n), n_slabs=n,
                           us_per_probability=1.0e6*t/len(Ls),
                           max_abs_error=float(np.max(np.abs(P - ref)))))
            print('  w=%.0e magnus  n=%-8d %.3e  %.1f us/prob'
                  % (width, n, mg[-1]['max_abs_error'],
                     mg[-1]['us_per_probability']), file=sys.stderr)

        case = {'width': width,
                'reference_unitarity': ref_unitarity,
                'reference_P_ee': [float(x) for x in ref[:, 0, 0]],
                'series': [{'name': 'Magnus', 'dial': 'n_slabs', 'points': mg}]}

        if HAVE_NPE:
            npe_pts = []
            for n in NPE_DIALS:
                P, used = npe_along_ray(width, n)
                t = timed(lambda n=n: npe_along_ray(width, n))
                npe_pts.append(dict(label=str(n), n_slabs=n, n_slabs_used=used,
                                    us_per_probability=1.0e6*t/len(Ls),
                                    max_abs_error=float(np.max(np.abs(P - ref)))))
                print('  w=%.0e npe     n=%-8d %.3e  %.1f us/prob'
                      % (width, n, npe_pts[-1]['max_abs_error'],
                         npe_pts[-1]['us_per_probability']), file=sys.stderr)
            case['series'].append({'name': 'NuOscProbExact', 'dial': 'n_slabs',
                                   'points': npe_pts})
            case['magnus_P_ee'] = [
                float(x) for x in magnus_along_ray(width, MG_DIALS[-1])[:, 0, 0]]
            case['npe_P_ee'] = [
                float(x) for x in npe_along_ray(width, NPE_DIALS[-1])[0][:, 0, 0]]
        out['cases'].append(case)
        checkpoint(out, 'case:%.0e' % width)

    for width in WIDTHS:
        if 'energy:%.0e' % width in done:
            continue
        out['prob_vs_energy'].append(prob_vs_energy(width))
        pe = out['prob_vs_energy'][-1]
        if pe['npe'] is not None:
            print('  w=%.0e prob-vs-E: worst |Magnus - NPE| = %.3e'
                  % (width, float(np.max(np.abs(np.array(pe['magnus'])
                                               - np.array(pe['npe']))))),
                  file=sys.stderr)
        checkpoint(out, 'energy:%.0e' % width)


if __name__ == '__main__':
    main()
