"""Issue #194: re-measure the timings quoted in the paper, on one quiet machine.

Every timed piece of code is read out of resources/paper/main.tex at run time -- a
listing by its \\label, a pysnip by a unique anchor string -- and run as printed, in a
namespace prepared the way the paper prepares it (the names an earlier listing defines).
Where the paper quotes a number but prints no code, the code here is built from the text
and the row says so.

    python docs/dev/measurements/issue194_paper_timings/run_all.py [--only 1,4,7] [--quick]

--quick runs every measurement once, with no warm-up, to check that the script works.
Prints one table (line | paper says | measured | suggested text) and writes the raw
timings to results_<hostname>_<date>.json next to this file.
"""

import os

# One thread per process, before NumPy is imported, as notebooks/gen_njobs_scaling.py
# does (and as the n_jobs figure of the paper states).  Subprocesses inherit it.
PINNED = ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
          'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
for _var in PINNED:
    os.environ[_var] = '1'
# The listings' warnings are not what is timed; joblib's workers do not inherit filters.
os.environ['PYTHONWARNINGS'] = 'ignore'

import argparse                                              # noqa: E402
import contextlib                                            # noqa: E402
import datetime                                              # noqa: E402
import io                                                    # noqa: E402
import json                                                  # noqa: E402
import math                                                  # noqa: E402
import pathlib                                               # noqa: E402
import platform                                              # noqa: E402
import socket                                                # noqa: E402
import statistics                                            # noqa: E402
import subprocess                                            # noqa: E402
import sys                                                   # noqa: E402
import time                                                  # noqa: E402
import warnings                                              # noqa: E402

import numpy as np                                           # noqa: E402

import magnus                                                # noqa: E402
import magnus.globaldefs as gd                               # noqa: E402
import magnus.hamiltonians as ham                            # noqa: E402
import magnus.oscprob as oscprob                             # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[3]
TEX_PATH = REPO/'resources'/'paper'/'main.tex'
TEX = TEX_PATH.read_text()

# Set in main(): repeats are multiplied down to 1 and warm-ups dropped under --quick.
QUICK = False
ROWS = []          # the table
RAW = {}           # every timing, for the JSON
CONTROL = []       # the stability control, (label, seconds)


# ------------------------------------------------------------------ reading the paper
def line_of(anchor):
    """The line number in main.tex of a string that must occur exactly once."""
    count = TEX.count(anchor)
    if count != 1:
        raise SystemExit('anchor found %d times in main.tex (needs 1): %r' % (count, anchor))
    return TEX[:TEX.index(anchor)].count('\n') + 1


def listing(label):
    """The code of the lstlisting carrying label={label}."""
    tag = 'label={%s}' % label
    line_of(tag)
    head_end = TEX.index('\n', TEX.index(tag)) + 1
    begin = TEX.rindex('\\begin{lstlisting}', 0, head_end)
    assert '\\end{lstlisting}' not in TEX[begin:head_end], label
    return TEX[head_end:TEX.index('\\end{lstlisting}', head_end)]


def pysnip(anchor):
    """The code of the pysnip block that contains `anchor`."""
    at = TEX.index(anchor)
    line_of(anchor)
    begin = TEX.rindex('\\begin{pysnip}', 0, at)
    end = TEX.index('\\end{pysnip}', at)
    assert '\\end{pysnip}' not in TEX[begin:at], anchor
    return TEX[TEX.index('\n', begin) + 1:end]


def cut(code, start=None, stop=None):
    """The part of `code` from the line holding `start` up to the line holding `stop`."""
    a = 0 if start is None else code.index(start)
    if start is not None:
        assert code.count(start) == 1, start
        a = code.rindex('\n', 0, a) + 1 if '\n' in code[:a] else 0
    if stop is None:
        return code[a:]
    assert code.count(stop) == 1, stop
    b = code.index(stop, a)
    b = code.rindex('\n', 0, b) + 1 if '\n' in code[a:b] else b
    return code[a:b]


def replace_once(code, old, new):
    assert code.count(old) == 1, 'expected exactly one %r' % old
    return code.replace(old, new)


def run_code(code, ns=None):
    """Execute paper code in a namespace (a fresh one unless given)."""
    ns = {} if ns is None else ns
    exec(compile(code, 'main.tex', 'exec'), ns)
    return ns


def as_function(code, ns):
    compiled = compile(code, 'main.tex', 'exec')
    return lambda: exec(compiled, ns)


# ------------------------------------------------------------------ timing
def reps(n):
    return 1 if QUICK else n


def measure(name, fn, n, warm=1, min_block=0.05):
    """Median of `n` timed blocks of `fn`, warm-up excluded.

    A block calls `fn` as many times as it takes to last `min_block` seconds (the
    paper's protocol, Sec. 7.1) and is divided by the number of calls; anything that
    takes longer than that is one call per block.  Returns the per-call seconds.
    """
    n = reps(n)
    for _ in range(0 if QUICK else warm):
        fn()
    calls = 1
    while True:
        t0 = time.perf_counter()
        for _ in range(calls):
            fn()
        elapsed = time.perf_counter() - t0
        if elapsed >= min_block or QUICK:
            break
        calls *= 2
    times = [elapsed/calls]
    while len(times) < n:
        t0 = time.perf_counter()
        for _ in range(calls):
            fn()
        times.append((time.perf_counter() - t0)/calls)
    RAW[name] = dict(seconds=times, calls_per_block=calls)
    print('    %-52s median %s  [%s-%s]  (%d x %d calls)'
          % (name, human(statistics.median(times)), human(min(times)), human(max(times)),
             len(times), calls), flush=True)
    return times


def silenced(fn):
    """`fn` with its stdout thrown away (for magnus's printed warnings)."""
    def quiet():
        with contextlib.redirect_stdout(io.StringIO()):
            return fn()
    return quiet


def sig_up(x):
    """`x` rounded up to one significant figure, for 'under ...'."""
    step = 10.0**math.floor(math.log10(x))
    return sig(math.ceil(x/step - 1e-9)*step, 1)


def human(s):
    return '%.3g ms' % (1e3*s) if s < 0.1 else '%.3g s' % s


def spread(times, unit='s'):
    k = 1e3 if unit == 'ms' else 1.0
    return '%.3g [%.3g-%.3g] %s' % (k*statistics.median(times), k*min(times),
                                    k*max(times), unit)


def sig(x, n):
    """`x` rounded to `n` significant figures, written without an exponent."""
    if x == 0 or not math.isfinite(x):
        return str(x)
    digits = n - 1 - math.floor(math.log10(abs(x)))
    x = round(x, digits)
    return '%.*f' % (max(digits, 0), x)


def row(anchor, paper, measured, suggested):
    ROWS.append((line_of(anchor), paper, measured, suggested))


def med(times):
    return statistics.median(times)


def control(label):
    """A fixed CPU-bound workload: best of 20 products of a 300x300 matrix."""
    a = np.random.default_rng(0).normal(size=(300, 300))
    best = math.inf
    for _ in range(20):
        t0 = time.perf_counter()
        a @ a
        best = min(best, time.perf_counter() - t0)
    CONTROL.append((label, best))
    print('  [control %-14s %.3f ms, drift %.3f]'
          % (label, 1e3*best, best/CONTROL[0][1]), flush=True)


# ------------------------------------------------------------------ 1. Listing 1
# A fresh process: starts Python, imports, then runs the listing read from stdin.
FRESH_CHILD = '''import time; t0 = time.perf_counter()
import sys, json, warnings
import numpy as np, magnus.globaldefs as gd, magnus.oscprob as oscprob
warnings.simplefilter('ignore')
t1 = time.perf_counter()
exec(compile(sys.stdin.read(), 'listing', 'exec'), {})
t2 = time.perf_counter()
print(json.dumps(dict(imports=t1 - t0, calls=t2 - t1)))
'''


def fresh(name, code, n):
    """Wall clock of a new process (with Python's start-up) and its calls alone."""
    runs = []
    for i in range((0 if QUICK else 1) + reps(n)):     # the first fills the kernel cache
        t0 = time.perf_counter()
        out = subprocess.run([sys.executable, '-c', FRESH_CHILD], input=code,
                             capture_output=True, text=True, check=True).stdout
        wall = time.perf_counter() - t0
        runs.append(dict(json.loads(out.strip().splitlines()[-1]), wall=wall))
    first, timed = (None, runs) if QUICK else (runs[0], runs[1:])
    RAW[name] = dict(first_run=first, runs=timed)
    walls = [r['wall'] for r in timed]
    calls = [r['calls'] for r in timed]
    print('    %-52s wall %s, calls %s' % (name, spread(walls), spread(calls)), flush=True)
    return walls, calls


def engines_of(code):
    """Which engine answered each of the four wrappers of Listing 1."""
    names = ['osc_prob_%dnu_matter_exp_density' % d for d in (2, 3, 4, 5)]
    real = {n: getattr(oscprob, n) for n in names}
    seen = []
    for n in names:
        def wrap(*a, _f=real[n], **k):
            info = {}
            out = _f(*a, strategy_info=info, **k)
            seen.append(info.get('engine'))
            return out
        setattr(oscprob, n, wrap)
    try:
        run_code(code)
    finally:
        for n in names:
            setattr(oscprob, n, real[n])
    return seen


def item1():
    printed = listing('lst:validation')
    default_tol = replace_once(printed, 'rtol=1e-12, atol=1e-14,', '')
    hybrid = replace_once(printed, 'magnus_exp_order=8)',
                          "magnus_exp_order=8, strategy='hybrid')")

    walls, calls = fresh('1 lst:validation, fresh process', printed, 5)
    d_walls, d_calls = fresh('1 lst:validation default tol, fresh process', default_tol, 5)
    warm = measure('1 lst:validation, warm', as_function(printed, {}), 5)
    d_warm = measure('1 lst:validation default tol, warm', as_function(default_tol, {}), 5)
    h_warm = measure('1 lst:validation strategy=hybrid, warm', as_function(hybrid, {}), 3)
    RAW['1 engines, as printed'] = engines_of(printed)
    RAW['1 engines, hybrid'] = engines_of(hybrid)
    print('    engines as printed %s; hybrid %s'
          % (RAW['1 engines, as printed'], RAW['1 engines, hybrid']))

    a = 'takes about half a second on a laptop'
    row(a, 'about half a second (warm)', spread(warm), 'about %s s' % sig(med(warm), 1))
    row(a, '(same, fresh process, calls only)', spread(calls),
        'about %s s' % sig(med(calls), 1))
    row(a, 'about a second incl. start of Python', spread(walls),
        'about %s s' % sig(med(walls), 1))
    a = 'the listing takes about 0.2~s, most of it'
    row(a, 'about 0.2 s at default tol (fresh, calls)', spread(d_calls),
        'about %s s' % sig(med(d_calls), 1))
    row(a, '(same, warm)', spread(d_warm), 'about %s s' % sig(med(d_warm), 1))
    row(a, '(same, fresh incl. start of Python)', spread(d_walls),
        'about %s s' % sig(med(d_walls), 1))
    row('some twenty times faster there than its adiabatic engine',
        'some twenty times faster than hybrid', 'hybrid %s; ratio %.1f'
        % (spread(h_warm), med(h_warm)/med(warm)),
        'some %s times' % sig(med(h_warm)/med(warm), 1))


# ------------------------------------------------------------------ 2. batched vs n_jobs
def njobs_setup():
    """The first pysnip of Sec. 'Running a scan in parallel': E, L, kw."""
    return pysnip('# batched engine takes this scan')


def item2():
    setup = cut(njobs_setup(), None, '# Every point has its own baseline')
    shared = pysnip('# engine takes it, at the default n_jobs = 1')
    shared = cut(shared, 'P = oscprob.osc_prob_3nu_earth(')
    ten = replace_once(shared, '**kw)', 'n_jobs=10, **kw)')
    results = {}
    for tol, sub in (('rtol=1e-6 (printed kw)', None),
                     ('rtol=atol=1e-3 (default)', ('rtol=1e-6, atol=1e-8', 'rtol=1e-3, atol=1e-3'))):
        s = setup if sub is None else replace_once(setup, *sub)
        ns = run_code(s)
        one = measure('2 batched, 5000 E, n_jobs=1, %s' % tol, as_function(shared, ns), 5)
        many = measure('2 per-point, 5000 E, n_jobs=10, %s' % tol, as_function(ten, ns), 3)
        for label, code in (('n_jobs=1', shared), ('n_jobs=10', ten)):
            ns['info'] = {}
            run_code(replace_once(code, '**kw)', 'strategy_info=info, **kw)'), ns)
            RAW['2 engine, %s, %s' % (label, tol)] = ns['info'].get('engine')
        print('    engines: n_jobs=1 %s, n_jobs=10 %s'
              % (RAW['2 engine, n_jobs=1, %s' % tol], RAW['2 engine, n_jobs=10, %s' % tol]))
        results[tol] = one, many
    a = 'the batched path takes 0.11~s in one process'
    for tol, (one, many) in results.items():
        row(a, '0.11 s in one process [%s]' % tol, spread(one), '%s~s' % sig(med(one), 2))
        row(a, 'ten processes take 1.1 s [%s]' % tol, spread(many),
            '%s~s' % sig(med(many), 2))


# ------------------------------------------------------------------ 3. batched vs per point
def item3():
    """Built from the text: 200-point scans, 2-5 flavors, batched against a loop.

    Energy scans use the paper's Es (0.1-10 GeV) and baseline scans its Ls (0-1300 km), at
    the L and E of lst:minimal; along the Earth, the chord at cos(theta_z) = -0.9 of the
    n_jobs section (energies at 11467.8 km; 200 baselines from 1 km to 11467.8 km at 1 GeV).
    """
    ns = run_code(listing('lst:minimal').split('P = oscprob')[0])      # E, L
    ns['np'] = np                       # the snippets use np; lst:minimal does not import it
    run_code(cut(pysnip('P = oscprob.osc_prob_3nu_vacuum(Es, L)'), None, 'P = '), ns)
    run_code(cut(pysnip('P = oscprob.osc_prob_3nu_vacuum(E, Ls)'), None, 'P = '), ns)
    E, L, Es, Ls = ns['E'], ns['L'], ns['Es'], ns['Ls']
    osc = gd.load_nufit_params('NuFIT 6.1')
    L_earth = 11467.8*gd.UNIT_KM
    Ls_earth = np.linspace(1.0, 11467.8, 200)*gd.UNIT_KM
    gains = {'vacuum': [], 'constant': [], 'earth': []}
    for d in (2, 3, 4, 5):
        # Sterile parameters left at the wrappers' defaults: with Listing 1's eV-scale
        # mixing the Earth energy scan is the paper's own exception (1.7x, same paragraph)
        par = dict(sth=osc['s13'], Dm2=osc['D31']) if d == 2 else dict(osc)
        cases = {
            'vacuum': (getattr(oscprob, 'osc_prob_%dnu_vacuum' % d), {},
                       [('E', Es, L), ('L', E, Ls)]),
            'constant': (getattr(oscprob, 'osc_prob_%dnu_matter_constant_density' % d),
                         dict(rho=3.0, density_matter_is_in_g_per_cm3=True),
                         [('E', Es, L), ('L', E, Ls)]),
            'earth': (getattr(oscprob, 'osc_prob_%dnu_earth' % d), dict(costhz=-0.9),
                      [('E', Es, L_earth), ('L', E, Ls_earth)])}
        for env, (f, extra, scans) in cases.items():
            kw = dict(par, **extra)
            for axis, e, l in scans:
                def batched(f=f, e=e, l=l, kw=kw):
                    return f(e, L=l, **kw)

                def looped(f=f, e=e, l=l, kw=kw):
                    es = np.broadcast_to(e, (200,))
                    ls = np.broadcast_to(l, (200,))
                    return [f(es[i], L=ls[i], **kw) for i in range(200)]
                tag = '3 %dnu %s, 200 %s' % (d, env, axis)
                tb = measure(tag + ', batched', batched, 5)
                tl = measure(tag + ', one call per point', looped, 3)
                gains[env].append(med(tl)/med(tb))
    RAW['3 gains'] = gains
    a = 'batching is 60--190 times faster'
    vc = gains['vacuum'] + gains['constant']
    row(a, '60-190x, vacuum and constant density',
        '%.0f-%.0fx' % (min(vc), max(vc)), '%s--%s times' % (sig(min(vc), 2), sig(max(vc), 2)))
    g = gains['earth']
    row(a, '11-120x, Earth chord', '%.0f-%.0fx' % (min(g), max(g)),
        '%s--%s times' % (sig(min(g), 2), sig(max(g), 2)))


# ------------------------------------------------------------------ 4. Listing lst:constant
def item4():
    code = listing('lst:constant')
    ns = run_code(code)                                   # every name the blocks use
    blocks = [
        ('wrappers, vacuum + matter, one energy', 'P_vac = oscprob.osc_prob_3nu_vacuum',
         '# Pme = 0.00377306', '0.090 ms', '# matter, in 0.090 ms'),
        ('wrapper, 200 energies', 'P_scan = oscprob.osc_prob_3nu_matter_constant_density(',
         '# P_scan.shape', '0.123 ms', 'in 0.123 ms'),
        ('scenario function, one energy', 'P_mat = oscprob.osc_prob_matter_std_potential(',
         '# The same number, to 1e-16', '31.9 ms', 'to 1e-16, in 31.9 ms'),
        ('osc_prob, vacuum + matter', 'P_vac = oscprob.osc_prob(H_vac/E, 0.0, L)',
         '# The same two numbers', '0.013 ms', 'The same two numbers, in 0.013 ms'),
        ('osc_prob_energy_baseline, 200 energies', 'P_scan = oscprob.osc_prob_energy_baseline(',
         '# The same 200 points', '3.6 ms', 'The same 200 points, in 3.6 ms')]
    t = {}
    for label, start, stop, paper, anchor in blocks:
        t[paper] = measure('4 lst:constant ' + label, as_function(cut(code, start, stop), ns), 7)
        n_sig = len(paper.split()[0].replace('0.', '').lstrip('0').replace('.', ''))
        row(anchor, '# ... in %s' % paper, spread(t[paper], 'ms'),
            'in %s ms' % sig(1e3*med(t[paper]), n_sig))
    a, b = med(t['0.123 ms']), med(t['0.090 ms'])
    row('Two hundred energies cost 0.123~ms vs.~0.090~ms', '0.123 ms vs. 0.090 ms',
        '%.3f vs %.3f ms' % (1e3*a, 1e3*b), '%s~ms vs.~%s~ms' % (sig(1e3*a, 3), sig(1e3*b, 2)))
    per = med(t['0.013 ms'])/2
    row('the cheapest route per probability, about 0.01~ms', 'about 0.01 ms per probability',
        '%.4f ms (block of two calls / 2)' % (1e3*per), 'about %s~ms' % sig(1e3*per, 1))
    e = med(t['3.6 ms'])
    row('For the 200 energies of Listing~\\ref{lst:constant}, that takes 3.6~ms',
        '3.6 ms', spread(t['3.6 ms'], 'ms'), '%s~ms' % sig(1e3*e, 2))
    row('in 0.123~ms, about thirty times faster', 'about thirty times faster',
        '%.1fx' % (e/a), 'about %s times' % sig(e/a, 1))


# ------------------------------------------------------------------ 5. Hamiltonians
def item5():
    ns = run_code(listing('lst:constant'))               # osc, E, Es, ham
    snip = pysnip('# Over 200 energies: 4.1 ms')
    run_code(cut(snip, None, '# Over 200 energies'), ns)  # H1, H0
    slow = measure('5 rebuild H at each of 200 energies',
                   as_function(cut(snip, 'Hs = [ham.hamiltonian_3nu_vacuum(e', '# ... and'), ns), 7)
    fast = measure('5 build once, divide by each energy',
                   as_function(cut(snip, 'Hs = [H0/e for e in Es]'), ns), 7)
    row('rebuilding the Hamiltonian at each energy takes about 4~ms', 'about 4 ms vs 0.2 ms',
        '%s vs %s' % (spread(slow, 'ms'), spread(fast, 'ms')),
        'about %s~ms ... %s~ms' % (sig(1e3*med(slow), 1), sig(1e3*med(fast), 1)))
    row('# Over 200 energies: 4.1 ms', '# 4.1 ms', spread(slow, 'ms'),
        '%s ms' % sig(1e3*med(slow), 2))
    row('# ... and 0.16 ms', '# 0.16 ms', spread(fast, 'ms'), '%s ms' % sig(1e3*med(fast), 2))

    # Above five flavors.  No code is printed, and the ratio depends strongly on the
    # energy (more slabs, more eigh), so it is measured at two: Listing 1's exponential
    # profile, one probability, nu_e -> nu_e, the default tolerance.  5 flavors through the
    # shipped builder with Listing 1's sterile mixing; 6 and 8 through h_vac_energy_indep,
    # the same matrix plus heavier states mixed with nu_e and nu_mu at s = 0.1.
    osc = gd.load_nufit_params('NuFIT 6.1')
    st5 = dict(s14=0.1**0.5, s24=0.1**0.5, s15=0.06**0.5, s25=0.06**0.5, D41=1.0, D51=1.7,
               s34=0.0, s35=0.0, d14=0.0, d15=0.0, d24=0.0, d35=0.0)
    h5 = ham.hamiltonian_5nu_vacuum_energy_independent(**osc, **st5)
    rho = lambda l: 3.e3*np.exp(-np.asarray(l)/(10.0*gd.UNIT_KM))   # noqa: E731
    common = dict(L0=0.0, nu_i=gd.NUE, nu_f=gd.NUE, density_matter_is_in_g_per_cm3=True)
    L5 = 25.0*gd.UNIT_KM

    def h_vac(n):
        """h5 with n - 5 heavier states (splittings 2.3, 2.9, ... where D51 = 1.7)."""
        w, v = np.linalg.eigh(h5)
        heavy = [w[-1]*(2.3 + 0.6*k)/1.7 for k in range(n - 5)]
        u = np.eye(n, dtype=complex)
        u[:5, :5] = v
        for k in range(5, n):
            for a in (0, 1):
                r = np.eye(n)
                r[a, a] = r[k, k] = math.sqrt(1 - 0.01)
                r[a, k], r[k, a] = 0.1, -0.1
                u = r @ u
        h = u @ np.diag(np.concatenate((w, heavy))) @ u.conj().T
        return 0.5*(h + h.conj().T)

    for e_gev in (0.05, 0.5):
        E5 = e_gev*gd.UNIT_GEV
        t = {5: measure('5 one probability, 5 flavors, %g GeV' % e_gev,
                        lambda: oscprob.osc_prob_matter_std_potential(
                            5, rho, E5, L5, dict(osc, **st5), **common), 7)}
        for n in (6, 8):
            hn = h_vac(n)
            # magnus prints two warnings per call above five flavors; kept off the screen
            t[n] = measure('5 one probability, %d flavors, %g GeV' % (n, e_gev),
                           silenced(lambda hn=hn, n=n: oscprob.osc_prob_matter_std_potential(
                               n, rho, E5, L5, osc, h_vac_energy_indep=hn, **common)), 7)
        r6, r8 = med(t[6])/med(t[5]), med(t[8])/med(t[5])
        row('three and a half times as long at eight',
            'twice at 6, 3.5x at 8 vs 5 [built here, %g GeV]' % e_gev,
            '%.2fx, %.2fx  (5 flavors: %s)' % (r6, r8, human(med(t[5]))),
            '%s ... %s times as long' % (sig(r6, 2), sig(r8, 2)))


# ------------------------------------------------------------------ 6. prob_vs, arrangement
ARRANGEMENT_RHO = {   # notebook 28's other three profiles; the listing holds the castle wall
    'serrated': np.tile(np.linspace(2.0, 8.0, 6), 4),
    'random wall': np.random.default_rng(20260801).permutation(
        np.where(np.arange(24) % 2 == 0, 2.0, 8.0)),
    'uniform': np.full(24, 5.0)}


def item6():
    code = listing('lst:prob_vs')
    ns = run_code(cut(code, None, '# Left column'))
    left = measure('6 lst:prob_vs left column (5000 L)',
                   as_function(cut(code, 'L = np.linspace(20.0, 500.0, 5000)*KM',
                                   '# Right column'), ns), 5)
    right = measure('6 lst:prob_vs right column (3000 E)',
                    as_function(cut(code, 'E = np.logspace(np.log10(3.0), 2.0, 3000)*MEV',
                                    '# The vacuum curves'), ns), 5)
    row('computes the whole left column in under half a second', 'left column: under 0.5 s',
        spread(left), 'in under %s s' % sig_up(med(left)))
    row('computes the whole right column in under a second', 'right column: under 1 s',
        spread(right), 'in under %s s' % sig_up(med(right)))

    code = listing('lst:arrangement')
    ns = run_code(cut(code, None, 'P, Pbar = ['))
    calls = as_function(cut(code, 'P, Pbar = ['), ns)
    castle = ns['rho'].copy()
    two = measure('6 lst:arrangement, castle wall, nu + nubar', calls, 5)

    def eight():
        for rho in [castle] + list(ARRANGEMENT_RHO.values()):
            ns['rho'] = rho
            calls()
        ns['rho'] = castle
    all8 = measure('6 lst:arrangement, four profiles, eight curves', eight, 5)
    row('take 0.1~s in all', 'eight curves take 0.1 s in all', spread(all8),
        '%s~s in all' % sig(med(all8), 1))
    row('take 0.1~s in all', '(the listing: two curves)', spread(two), '')


# ------------------------------------------------------------------ 7. the Sun
def item7():
    code = listing('lst:sun')
    ns = run_code(cut(code, None, '# Standard three-flavor oscillations'))
    parts = [('3 flavors', 'P3 = oscprob.osc_prob_3nu_sun', '# The same, with non-standard'),
             ('3 flavors + NSI', 'eps = dict(', '# One and two sterile states'),
             ('3+1', 'ster = dict(', 'P5 = oscprob'),
             ('3+2', 'P5 = oscprob', '# P3 = 0.5449')]
    t = []
    for label, start, stop in parts:
        t.append(measure('7 lst:sun ' + label, as_function(cut(code, start, stop), ns), 3))
    meds = [med(x) for x in t]
    row('take about 0.5~s, 0.5~s, 0.8~s, and 1.2~s', '0.5, 0.5, 0.8, 1.2 s',
        ', '.join('%.2f' % m for m in meds) + ' s',
        ', '.join('%s~s' % sig(m, 1 if m < 1 else 2) for m in meds))
    row('The four curves take about three seconds', 'about three seconds',
        '%.2f s (sum of medians)' % sum(meds), 'about %s seconds' % sig(sum(meds), 1))

    snip = pysnip('# P = 0.666, 0.092, 0.259 in about 7 s')
    slow = measure('7 instantaneous solar snippet', as_function(snip, ns), 3)
    row('# P = 0.666, 0.092, 0.259 in about 7 s', 'in about 7 s', spread(slow),
        'in about %s s' % sig(med(slow), 1))


# ------------------------------------------------------------------ 8. n_jobs speed-up
def item8():
    snip = njobs_setup()
    out = {}
    for n_points, n in ((500, 5), (2000, 3), (5000, 3)):
        code = replace_once(snip, 'N = 5000', 'N = %d' % n_points)
        ns = run_code(cut(code, None, 'P = oscprob'))
        call = cut(code, 'P = oscprob.osc_prob_3nu_earth(')
        fns = {k: as_function(replace_once(call, 'n_jobs=4', 'n_jobs=%d' % k), ns) for k in (1, 4)}
        times = {1: [], 4: []}
        for k in (1, 4):                    # warm-up: kernels, and joblib's worker pool,
            if not QUICK:                   # which n_jobs=1 leaves alone between rounds
                fns[k]()
        for _ in range(reps(n)):            # interleaved rounds
            for k in (1, 4):
                t0 = time.perf_counter()
                fns[k]()
                times[k].append(time.perf_counter() - t0)
        RAW['8 n_jobs, %d points' % n_points] = times
        out[n_points] = med(times[1])/med(times[4])
        print('    8 %5d points: n_jobs=1 %s, n_jobs=4 %s, speed-up %.2f'
              % (n_points, spread(times[1]), spread(times[4]), out[n_points]), flush=True)
    row('four workers return $2.4$ times the speed', '2.4 at 5000, 2.3 at 2000, 1.8 at 500',
        '%.2f, %.2f, %.2f' % (out[5000], out[2000], out[500]),
        '$%.1f$ times ... $%.1f$ at $2\\,000$ and $%.1f$ at $500$'
        % (out[5000], out[2000], out[500]))


# ------------------------------------------------------------------ 9. the stack
def item9():
    def read(path, key):
        try:
            for ln in open(path):
                if ln.startswith(key):
                    return ln.split(':', 1)[1].strip() if ':' in ln else ln.split('=', 1)[1].strip('"\n')
        except OSError:
            pass
        return '?'
    import numba
    import scipy
    import joblib
    stack = dict(python=platform.python_version(), numpy=np.__version__,
                 numba=numba.__version__, scipy=scipy.__version__, joblib=joblib.__version__,
                 magnus=magnus.__version__, kernel=platform.release(),
                 os=read('/etc/os-release', 'PRETTY_NAME'),
                 cpu=read('/proc/cpuinfo', 'model name'), logical_cpus=os.cpu_count(),
                 memory=read('/proc/meminfo', 'MemTotal'))
    RAW['9 stack'] = stack
    row('running Ubuntu 24.04 LTS on Linux kernel', 'i5-1334U, Ubuntu 24.04, kernel 7.0.0, '
        'Python 3.12.7, NumPy 1.26.4, magnus 1.1.1',
        '; '.join('%s %s' % kv for kv in stack.items()),
        'running %s on Linux kernel %s, with Python %s, {\\tt NumPy} %s, and \\magnus\\ %s'
        % (stack['os'], stack['kernel'].split('-')[0], stack['python'], stack['numpy'],
           stack['magnus']))


ITEMS = {1: item1, 2: item2, 3: item3, 4: item4, 5: item5, 6: item6, 7: item7, 8: item8,
         9: item9}


def main():
    global QUICK
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    parser.add_argument('--only', help='comma-separated item numbers, 1-9')
    parser.add_argument('--quick', action='store_true',
                        help='one run of everything, no warm-up: checks the script only')
    args = parser.parse_args()
    QUICK = args.quick
    chosen = [int(x) for x in args.only.split(',')] if args.only else sorted(ITEMS)
    warnings.simplefilter('ignore')        # the listings' own warnings are not timed apart

    print('issue #194 paper timings%s; threads pinned: %s=1; load average %.2f'
          % (' (QUICK: numbers mean nothing)' if QUICK else '', '/'.join(PINNED),
             os.getloadavg()[0]), flush=True)
    started = time.perf_counter()
    control('start')
    for k in chosen:
        print('item %d' % k, flush=True)
        t0 = time.perf_counter()
        ITEMS[k]()
        RAW['item %d seconds' % k] = time.perf_counter() - t0
        control('after item %d' % k)

    drift = [c/CONTROL[0][1] for _, c in CONTROL]
    print('\n| line | paper says | measured (median [min-max]) | suggested text |')
    print('|---|---|---|---|')
    for r in ROWS:
        print('| %d | %s | %s | %s |' % r)
    print('\nstability control (300x300 matmul, best of 20): start %.3f ms; drift ratio '
          'min %.3f, max %.3f%s'
          % (1e3*CONTROL[0][1], min(drift), max(drift),
             '  <-- above 1.10: the machine was busy, rerun' if max(drift) > 1.10 else ''))
    print('threads: %s set to 1 for this process and its children; load average now %.2f; '
          'total %.0f s' % (', '.join(PINNED), os.getloadavg()[0],
                            time.perf_counter() - started))

    out = HERE/('results_%s_%s.json' % (socket.gethostname(), datetime.date.today().isoformat()))
    json.dump(dict(quick=QUICK, items=chosen, pinned={v: '1' for v in PINNED},
                   tex=str(TEX_PATH.relative_to(REPO)), control=CONTROL, rows=ROWS, raw=RAW,
                   platform=platform.platform()),
              open(out, 'w'), indent=1, default=str)
    print('raw timings written to %s' % out)

    print('\nIssue #194 part 2 (notebook 28 timings cache), only if wanted:\n'
          '  pip install -e ".[notebooks]"\n'
          '  env -u MAGNUS_PAPER_CACHE_ONLY MAGNUS_PAPER_RETIME=1 '
          'python notebooks/make_notebooks.py --only 28\n'
          '  (rewrites notebooks/28_magnus_paper_figures.ipynb and '
          'notebooks/paper_figure_cache.json; commit both)')


if __name__ == '__main__':
    main()
