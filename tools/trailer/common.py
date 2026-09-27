r"""Shared pieces of the trailer's generation scripts (issue #99).

Paths, the look (colors and fonts), the two physics setups the new scenes are built on, and an
independent reference integrator used to check what Magnus returns for them.  Nothing here edits
a figure source: every plot the trailer shows is remade by the scripts in this directory, and
everything they write goes to ``tools/trailer/build/``, which git ignores.
"""
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parents[1]
BUILD = HERE / 'build'
sys.path.insert(0, str(REPO / 'src'))

import magnus.globaldefs as gd  # noqa: E402
import magnus.matter as matter  # noqa: E402

FPS, BPM = 30, 96
SPB = 60.0 / BPM                        # seconds per beat

# The look: one dark theme, three accents (nu_e blue, nu_mu amber, nu_tau teal) and two more.
BG, PANEL, LINE, INK, MUT = '#06080d', '#0f1626', '#232c3f', '#e9edf5', '#a3acbf'   # MUT: labels and ticks
BLUE, AMBER, TEAL, ROSE, VIOLET = '#5aa9ff', '#f0a33c', '#3fd0a4', '#ff6b8b', '#b48cff'
# One typeface throughout, titles to tick labels to mathematics: Inter, a lean sans-serif with Greek.
# MONO and DISP are kept as names for the two roles, and both are Inter.
MONO = DISP = 'Inter'

# Inter is under the SIL Open Font License.  Its static files come from the official release,
# downloaded once on first use; Matplotlib renders only the default instance of a variable font,
# so the static weights are needed for bold and italic.
FONTS = HERE / 'fonts'                   # Inter 4.1, from github.com/rsms/inter/releases (extras/ttf)
INTER_FILES = ('Regular', 'Italic', 'Medium', 'SemiBold', 'Bold', 'BoldItalic')


def setup_matplotlib():
    """Agg backend, Inter registered and used for everything, the mathematics included.  The font
    files ship with the trailer, in ``fonts/`` (Inter 4.1, SIL Open Font License)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    wanted = [FONTS / ('Inter-%s.ttf' % w) for w in INTER_FILES]
    for f in wanted:
        fm.fontManager.addfont(str(f))
    plt.rcParams.update({'font.family': 'Inter', 'mathtext.fontset': 'custom', 'mathtext.rm': 'Inter',
                         'mathtext.it': 'Inter:italic', 'mathtext.bf': 'Inter:bold', 'mathtext.sf': 'Inter', 'mathtext.cal': 'Inter',
                         'axes.unicode_minus': True})
    return plt


def ease(x):
    """Smoothstep on [0, 1]."""
    x = float(np.clip(x, 0.0, 1.0))
    return x * x * (3 - 2 * x)


def seg(u, a, b):
    """Progress of a sub-animation that runs from u = a to u = b."""
    return ease((u - a) / (b - a))


# ------------------------------------------------------------------ the opening's two tracks
OPENING_E_GEV = 3.0                      # near enough to the MSW resonance that matter shows
OPENING_L_KM = np.linspace(0.0, 10000.0, 1001)


def opening_rho(l):
    """Density along the opening's matter track [g/cm^3], ``l`` in natural units: a smooth rise
    and fall from 3 to 10 g/cm^3 with a ripple.  Illustrative, not a real profile."""
    x = np.asarray(l, float) / gd.UNIT_KM / 10000.0
    return 3.0 + 7.0 * np.exp(-((x - 0.5) / 0.18)**2) + 0.6 * np.sin(2 * np.pi * 6 * x)


# ------------------------------------------------------------------ the adiabatic scene
# Two flavors at 10 MeV.  Lengths in units of 1/DELTA, DELTA = Dm2/2E (1/DELTA = 52.6 km);
# the potential in units of DELTA.  The envelope crosses the MSW resonance slowly, a tanh shock
# of width W crosses it fast, and a second slow decline crosses it again.  S stretches the
# envelope, not the shock.  S = 30, W = 2 is the scene: strategy='auto' answers long paths with
# the hybrid, which finds the shock window itself (issue #100 has a case where it cannot).
ADIAB_E, ADIAB_TH, ADIAB_DM2 = 10.e6, 0.3, 7.5e-5          # eV, rad, eV^2
ADIAB_S, ADIAB_W = 30.0, 2.0
DELTA = ADIAB_DM2 / (2 * ADIAB_E)
C2, S2 = np.cos(2 * ADIAB_TH), np.sin(2 * ADIAB_TH)


def v_over_delta(x, S=ADIAB_S, W=ADIAB_W):
    """Matter potential over DELTA; x in units of 1/DELTA."""
    x = np.asarray(x, float)
    xs = 1500.0 * S
    up = 0.5 * (1 + np.tanh((x - xs) / W))
    return 3.0 * np.exp(-x / (500.0 * S)) + 2.0 * up * np.exp(-(x - xs).clip(0) / (400.0 * S))


def adiab_rho_func(S=ADIAB_S, W=ADIAB_W):
    """The same profile as a density [g/cm^3] of natural-unit position, for the Magnus wrappers."""
    rho1 = DELTA / matter.vcc_func_from_rho_func(1.0, density_matter_is_in_g_per_cm3=True)
    return lambda l: rho1 * v_over_delta(np.asarray(l) * DELTA, S, W)


def reference_pee(L, S=ADIAB_S, W=ADIAB_W, n=2**21):
    """P(nu_e -> nu_e) over [0, L] (units of 1/DELTA), independent of Magnus: the exponential
    midpoint rule with closed-form 2x2 steps, dense around the shock, multiplied by a pairwise
    tree.  Second order; at n = 2^21 it agrees with 2^23 to seven digits on this scene."""
    xs = 1500.0 * S
    fine = np.clip(np.linspace(xs - 30 * W, xs + 30 * W, n // 8 + 1), 0, L)
    e = np.unique(np.concatenate([np.linspace(0, L, n + 1), fine]))
    m, h = 0.5 * (e[1:] + e[:-1]), np.diff(e)
    v = v_over_delta(m, S, W)
    a, bx, bz = 0.5 * v, 0.5 * S2 * np.ones_like(v), 0.5 * (v - C2)
    nb = np.hypot(bx, bz)
    c, s, g = np.cos(nb * h), np.sin(nb * h) / nb, np.exp(-1j * a * h)
    U = np.empty((len(h), 2, 2), complex)
    U[:, 0, 0], U[:, 1, 1] = g * (c - 1j * s * bz), g * (c + 1j * s * bz)
    U[:, 0, 1] = U[:, 1, 0] = -1j * g * s * bx
    while len(U) > 1:
        if len(U) % 2:
            U = np.concatenate([U, np.eye(2)[None]])
        U = U[1::2] @ U[0::2]
    return abs(U[0][0, 0])**2
