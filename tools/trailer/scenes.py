r"""The trailer's new scenes, each drawn as a function of its progress u in [0, 1].

A still is the scene at u = 1; an animation is the scene sampled at 30 fps over the shot's
length (see ``render.py``).  The numbers come from ``build/`` (see ``data.py``); nothing here
computes physics.
"""
import json

import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

from common import (BUILD, BG, PANEL, LINE, INK, MUT, BLUE, AMBER, TEAL, ROSE, VIOLET, MONO, DISP, seg)

_cache = {}


def load(name):
    """build/<name>.npz or .json, read once."""
    if name not in _cache:
        path = BUILD / name
        _cache[name] = json.loads(path.read_text()) if name.endswith('.json') else dict(np.load(path))
    return _cache[name]


def box(ax, x, y, w, h, fc=PANEL, ec=LINE, lw=1.2, r=0.12, alpha=1.0):
    if alpha > 0:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=%g' % r,
                                    fc=fc, ec=ec, lw=lw, alpha=alpha))


def arrow(ax, a, b, col=MUT, lw=1.6, alpha=1.0, p=1.0):
    """An arrow from a toward b, drawn to fraction p of its length."""
    if alpha > 0 and p > 0.02:
        a = np.asarray(a, float)
        b = a + p * (np.asarray(b, float) - a)
        ax.add_patch(FancyArrowPatch(tuple(a), tuple(b), arrowstyle='-|>', mutation_scale=16, color=col,
                                     lw=lw, alpha=alpha))


def t(ax, x, y, s, size=14, col=INK, fam=MONO, alpha=1.0, **kw):
    if alpha > 0:
        ax.text(x, y, s, fontsize=size, color=col, fontfamily=fam, alpha=alpha, **kw)


def plain(ax):
    ax.set_facecolor(BG)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_yticks([])
    ax.tick_params(colors=MUT, labelsize=11, length=0)


def monoticks(ax):
    for tk in ax.get_xticklabels() + ax.get_yticklabels():
        tk.set_fontfamily(MONO)


def glow(ax, x, y, col, sizes=((1400, 0.10), (600, 0.22), (240, 0.45), (90, 1.0))):
    for s, a in sizes:
        ax.scatter([x], [y], s=s, color='#ffffff' if a == 1.0 else col, alpha=a, lw=0, zorder=6, clip_on=False)


# ------------------------------------------------------------------ act 1: two journeys
def opening(fig, ax, u):
    """One neutrino rides two tracks, vacuum above and matter below; each track's P(nu_mu -> nu_e)
    draws itself up to where the neutrino is."""
    d = load('opening.npz')
    L, rho, Pv, Pm = d['L'], d['rho'], d['Pv'], d['Pm']
    k = int(round((0.02 + 0.98 * np.clip(u / 0.92, 0, 1)) * (len(L) - 1)))
    top = 1.12 * max(Pv.max(), Pm.max())
    cmap = LinearSegmentedColormap.from_list('rho', ['#1a1410', '#6b3d12', AMBER])
    for y0, name, P, col, dens in ((0.56, 'VACUUM', Pv, BLUE, None), (0.10, 'MATTER', Pm, AMBER, rho)):
        a = fig.add_axes([0.08, y0, 0.86, 0.30])
        plain(a)
        a.set_xlim(0, L[-1])
        a.set_ylim(-0.26 * top, top * 1.1)
        if dens is None:
            a.axhspan(-0.24 * top, -0.08 * top, color='#101a2c', lw=0)
        else:
            a.imshow(dens[None, :], extent=(0, L[-1], -0.24 * top, -0.08 * top), aspect='auto', cmap=cmap,
                     vmin=dens.min(), vmax=dens.max(), zorder=1)
        a.plot(L, P, color=col, lw=1.2, alpha=0.13)
        a.plot(L[:k + 1], P[:k + 1], color=col, lw=3.2, solid_capstyle='round')
        a.fill_between(L[:k + 1], 0, P[:k + 1], color=col, alpha=0.10, lw=0)
        glow(a, L[k], -0.16 * top, col)
        a.scatter([L[k]], [P[k]], s=60, color=col, zorder=7, lw=0)
        a.text(0, top * 1.08, name, color=col, fontfamily=DISP, fontsize=17, weight=700, va='top')
        a.text(L[-1], top * 1.08, r'$P(\nu_\mu \to \nu_e)$ = %.2f' % P[k], color=INK, fontfamily=MONO,
               fontsize=17, ha='right', va='top')
        a.set_xticks([0, 2500, 5000, 7500, 10000] if dens is not None else [])
        monoticks(a)
        if dens is not None:
            a.set_xlabel('distance travelled  L  [km]', color=MUT, fontfamily=MONO, fontsize=13)
    t(ax, 15.04, 0.4, 'E = %g GeV' % d['E'], 13, MUT, ha='right')


# ------------------------------------------------------------------ act 4: Magnus decides
def switch(fig, ax, u):
    """The neutrino crosses the MSW resonance three times; the levels draw along its path, the
    slow crossings in teal (adiabatic transport), the shock crossing amber (Magnus patch, found by
    Magnus and magnified at right); P_ee lands one dot per call."""
    d = load('adiabatic.npz')
    x, lam, Ls, P, res, wins, c2, km = (d[k] for k in ('x', 'lam', 'Ls', 'P', 'res', 'wins', 'c2', 'km'))
    engines = d['engines']
    X, LX = x * km / 1e6, Ls * km / 1e6
    xn = X[-1] * np.clip(u / 0.9, 0.004, 1)
    wa, wb = wins[0] * km / 1e6
    a = fig.add_axes([0.07, 0.40, 0.60, 0.44])
    plain(a)
    a.set_xlim(0, X[-1])
    a.set_xticks([])
    k = np.searchsorted(X, xn)
    for j in range(2):
        a.plot(X, lam[:, j], color=TEAL, lw=1.2, alpha=0.12)
        a.plot(X[:k + 1], lam[:k + 1, j], color=TEAL, lw=2.4)
    lo, hi = lam.min() - 0.05, lam.max() + 0.05
    a.set_ylim(lo, hi)
    hit = seg(xn, wa, wa + 0.12)                      # the window lights once the neutrino reaches it
    if hit > 0:
        a.axvspan(wa - 0.02, wb + 0.02, color=AMBER, alpha=0.9 * hit, lw=0)
    for r in res:
        if r * km / 1e6 <= xn:
            a.axvline(r * km / 1e6, color=MUT, lw=0.8, ls=(0, (2, 4)))
    glow(a, xn, lam[min(k, len(X) - 1), 1], BLUE)
    t(a, 0, hi, 'THE LEVELS THE NEUTRINO FOLLOWS', 14, TEAL, DISP, weight=700, va='top')
    labels = ((res[0], 'slow crossing\nadiabatic transport', TEAL, 'left', lo),
              (wins[0][0], 'shock crossing  \nnon-adiabatic  ', AMBER, 'right', hi * 0.93 + lo * 0.07),
              (res[2], 'slow crossing\nadiabatic transport', TEAL, 'left', lo))
    for pos, s, col, ha, y in labels:
        p = seg(xn, pos * km / 1e6, pos * km / 1e6 + 0.15)
        t(a, pos * km / 1e6 + (0.03 if ha == 'left' else 0), y, s, 12, col, ha=ha,
          va='bottom' if y == lo else 'top', alpha=p)
    # the magnifier opens as the neutrino reaches the shock
    m = seg(xn, wa, wa + 0.35)
    if m > 0:
        ai = fig.add_axes([0.71, 0.40, 0.25, 0.44])
        plain(ai)
        S_, W_ = float(d['S']), float(d['W'])
        xf = np.linspace(wins[0][0] - 25, wins[0][1] + 25, 2001)
        xs, s2 = 1500.0 * S_, float(d['s2'])
        vf = 3.0 * np.exp(-xf / (500 * S_)) + 2.0 * 0.5 * (1 + np.tanh((xf - xs) / W_)) * np.exp(-(xf - xs).clip(0) / (400 * S_))
        lamf = np.stack([vf / 2 - np.hypot((vf - c2) / 2, s2 / 2), vf / 2 + np.hypot((vf - c2) / 2, s2 / 2)], 1)
        xx = (xf - xf[0]) * km
        ai.axvspan((wins[0][0] - xf[0]) * km, (wins[0][1] - xf[0]) * km, color=AMBER, alpha=0.18 * m, lw=0)
        for j in range(2):
            ai.plot(xx, lamf[:, j], color=AMBER, lw=2.4, alpha=m)
        ai.set_xlim(xx[0], xx[-1])
        ai.set_xticks([])
        y0 = ai.get_ylim()[0]
        ai.plot([xx[-1] - 1100, xx[-1] - 100], [y0, y0], color=MUT, lw=2, clip_on=False, alpha=m)
        ai.text(xx[-1] - 600, y0, '1000 km\n', color=MUT, fontfamily=MONO, fontsize=11, ha='center', va='bottom', alpha=m)
        ai.text(0.03, 0.97, 'MAGNUS PATCH', transform=ai.transAxes, color=AMBER, fontfamily=DISP, fontsize=14,
                weight=700, va='top', alpha=m)
        ai.text(0.03, 0.87, 'window found by Magnus,\nsolved exactly', transform=ai.transAxes, color=INK,
                fontfamily=MONO, fontsize=12, va='top', alpha=m)
        for s in ai.spines.values():
            s.set_visible(True)
            s.set_color(AMBER)
            s.set_alpha(m)
    # one dot per call, the engine that answered under the dot's path length
    ap = fig.add_axes([0.07, 0.09, 0.60, 0.24])
    plain(ap)
    ap.set_xlim(0, X[-1])
    ap.set_ylim(-0.05, 1.3)
    kk = np.searchsorted(LX, xn)
    ap.scatter(LX[kk:], P[kk:], s=10, color=BLUE, alpha=0.12, lw=0)
    if kk:
        ap.scatter(LX[:kk], P[:kk], s=22, color=BLUE, lw=0)
        j = kk - 1
        ap.scatter([LX[j]], [P[j]], s=90, color='#fff', zorder=5, lw=0)
        ap.text(0, 1.3, r'$P(\nu_e \to \nu_e)$ = %.2f     one dot, one call:  %s' % (P[j], engines[j]), color=INK,
                fontfamily=MONO, fontsize=15, va='top')
    ap.set_xticks([0, 1, 2, 3, 4])
    ap.set_xticklabels(['0', '1', '2', '3', '4 million km'])
    monoticks(ap)


# ------------------------------------------------------------------ diagram 1: the burden
def burden(fig, ax, u):
    rows = ['Reactor', 'Accelerator', 'Through the Earth', 'The Sun', 'Supernova']
    cols = ['Standard', 'Sterile', 'NSI', 'Lorentz viol.', 'Long-range']
    x0, y0, cw, ch = 3.1, 1.2, 1.55, 1.12
    t(ax, x0, 8.35, 'EVERY SETUP × EVERY THEORY', 15, AMBER, DISP, weight=700, alpha=seg(u, 0, 0.08))
    for j, c in enumerate(cols):
        t(ax, x0 + j * cw + cw / 2, y0 + 5 * ch + 0.25, c, 12.5, MUT, ha='center', alpha=seg(u, 0.02 + 0.02 * j, 0.1 + 0.02 * j))
    fade = 1 - 0.55 * seg(u, 0.62, 0.8)                 # the grid dims as it collapses
    for i, r in enumerate(rows):
        t(ax, x0 - 0.2, y0 + (4 - i) * ch + ch / 2, r, 12.5, MUT, ha='right', va='center',
          alpha=seg(u, 0.02 + 0.02 * i, 0.1 + 0.02 * i))
        for j in range(5):
            n = i * 5 + j
            a = seg(u, 0.08 + 0.02 * n, 0.12 + 0.02 * n) * fade
            xx, yy = x0 + j * cw + 0.1, y0 + (4 - i) * ch + 0.12
            box(ax, xx, yy, cw - 0.2, ch - 0.24, fc='#141b2b', ec='#2d3850', r=0.06, alpha=0.9 * a)
            t(ax, xx + (cw - 0.2) / 2, yy + (ch - 0.24) / 2, 'solver_%02d.py' % (n + 1), 9.5, AMBER, ha='center',
              va='center', alpha=0.85 * a)
    for i in range(5):
        arrow(ax, (x0 + 5 * cw + 0.05, y0 + i * ch + ch / 2), (12.25, 4.0), col=AMBER, lw=1.0, alpha=0.5,
              p=seg(u, 0.62 + 0.02 * i, 0.76 + 0.02 * i))
    b = seg(u, 0.74, 0.86)
    box(ax, 12.3, 3.1, 3.1, 1.8, fc='#1a1406', ec=AMBER, lw=2.2, r=0.2, alpha=b)
    t(ax, 13.85, 4.25, r'$H(E,\,x)$', 30, INK, ha='center', va='center', alpha=b)
    t(ax, 13.85, 3.45, 'one solver', 13, AMBER, ha='center', va='center', alpha=seg(u, 0.84, 0.92))
    t(ax, 13.85, 2.5, '25 files → 1 function', 12, MUT, ha='center', alpha=seg(u, 0.9, 0.98))


# ------------------------------------------------------------------ diagram 2: the Hamiltonian
def hamiltonian(fig, ax, u):
    D = load('diagrams.json')
    terms = [(r'$H(E,\,x) \;=$', 3.2, INK, 0.0), (r'$\dfrac{H_{\rm vac}}{E}$', 6.05, BLUE, 0.12),
             (r'$+\; V(x)\,P_e$', 8.6, AMBER, 0.38), (r'$+\; H_{\rm new}(x)$', 12.0, VIOLET, 0.64)]
    for s, x, col, a in terms:
        t(ax, x, 7.7, s, 34, col, ha='center', va='center', alpha=seg(u, a, a + 0.1))
    cards = [(1.0, BLUE, 'VACUUM', 'mixing + mass splittings', 0.16), (6.0, AMBER, 'MATTER', 'constant or not', 0.42),
             (11.0, VIOLET, 'BEYOND THE SM', 'your new physics', 0.68)]
    for x, col, head, sub, a in cards:
        v = seg(u, a, a + 0.1)
        box(ax, x, 1.1, 4.2, 5.0, ec=col, lw=1.8, alpha=v)
        t(ax, x + 0.3, 5.6, head, 15, col, DISP, weight=700, alpha=v)
        t(ax, x + 0.3, 5.1, sub, 12.5, MUT, alpha=v)
    U3 = np.asarray(D['U']['3'])
    for k, yl in enumerate((2.0, 2.5, 4.0)):
        g, left = seg(u, 0.2 + 0.04 * k, 0.34 + 0.04 * k), 1.6
        for f, col in enumerate((BLUE, AMBER, TEAL)):
            w = 3.0 * U3[f, k] * g
            if w > 0:
                ax.add_patch(Rectangle((left, yl), w, 0.28, fc=col, ec='none'))
            left += w
        t(ax, 1.45, yl + 0.14, r'$\nu_%d$' % (k + 1), 12, MUT, ha='right', va='center', alpha=g)
    t(ax, 1.6, 1.45, r'$\nu_e$  $\nu_\mu$  $\nu_\tau$  content of each mass state', 10.5, MUT, alpha=seg(u, 0.3, 0.38))
    xs = np.linspace(6.4, 9.8, 200)
    g = seg(u, 0.46, 0.62)
    k = max(2, int(g * len(xs)))
    if g > 0:
        ax.plot(xs[:k], 2.3 + 0 * xs[:k], color=AMBER, lw=2.4, alpha=0.55)
        ax.plot(xs[:k], (2.9 + 1.2 * np.exp(-((xs - 8.1) / 0.7)**2) + 0.12 * np.sin(9 * xs))[:k], color=AMBER, lw=2.6)
    t(ax, 9.9, 2.3, 'constant', 11, MUT, va='center', alpha=seg(u, 0.58, 0.64))
    t(ax, 9.3, 4.35, 'varying', 11, MUT, va='center', alpha=seg(u, 0.58, 0.64))
    for k, c in enumerate(['sterile neutrinos', 'non-standard interactions', 'Lorentz violation', 'long-range forces']):
        v = seg(u, 0.72 + 0.05 * k, 0.8 + 0.05 * k)
        box(ax, 11.4, 3.9 - 0.62 * k, 3.4, 0.46, fc='#161028', ec=VIOLET, lw=1.2, r=0.23, alpha=v)
        t(ax, 13.1, 4.13 - 0.62 * k, c, 11.5, VIOLET, ha='center', va='center', alpha=v)
    t(ax, 8, 0.45, 'any sum of these; any number of flavors', 13, MUT, ha='center', alpha=seg(u, 0.92, 1.0))


# ------------------------------------------------------------------ diagram 3: how Magnus works
def slabs(fig, ax, u):
    D = load('diagrams.json')
    d = load('opening.npz')
    t(ax, 0.8, 8.2, 'SLICE THE PATH. MULTIPLY. REFINE UNTIL TWO ANSWERS AGREE.', 16, TEAL, DISP, weight=700,
      alpha=seg(u, 0, 0.08))
    X = 0.8 + 9.4 * d['L'] / d['L'][-1]
    Y = 5.0 + 1.9 * (d['rho'] - d['rho'].min()) / (d['rho'].max() - d['rho'].min())
    k = max(2, int(seg(u, 0.0, 0.2) * len(X)))
    ax.fill_between(X[:k], 5.0, Y[:k], color=AMBER, alpha=0.12, lw=0)
    ax.plot(X[:k], Y[:k], color=AMBER, lw=2)
    t(ax, 0.8, 7.35, 'density along the path', 11, MUT, alpha=seg(u, 0.05, 0.15))
    for j in range(7):
        v = seg(u, 0.18 + 0.02 * j, 0.24 + 0.02 * j)
        if v > 0:
            ax.plot([0.8 + 9.4 * j / 6] * 2, [4.85, 4.85 + 2.2 * v], color=TEAL, lw=1.3)
    for j in range(6):
        t(ax, 0.8 + 9.4 * (j + 0.5) / 6, 4.45, r'$e^{\Omega_%d}$' % (j + 1), 16, TEAL, ha='center', va='center',
          alpha=seg(u, 0.3 + 0.02 * j, 0.36 + 0.02 * j))
    t(ax, 5.5, 3.55, r'$U \;=\; e^{\Omega_6}\,e^{\Omega_5}\cdots e^{\Omega_1}$', 22, INK, ha='center', va='center',
      alpha=seg(u, 0.44, 0.52))
    t(ax, 5.5, 2.85, r'$P = |U_{e\mu}|^2$, unitary by construction', 13, MUT, ha='center', alpha=seg(u, 0.5, 0.56))
    lad = D['ladder']
    t(ax, 11.0, 7.35, 'slabs    P(νμ→νe)', 13, MUT, alpha=seg(u, 0.5, 0.55))
    for k, (ns, P) in enumerate(lad):
        v = seg(u, 0.54 + 0.055 * k, 0.58 + 0.055 * k)
        yk = 6.75 - 0.62 * k
        ok = k >= 1 and abs(P - lad[k - 1][1]) <= 1e-5 + 1e-5 * abs(lad[k - 1][1])
        t(ax, 11.0, yk, '%5d    %.5f' % (ns, P), 15, TEAL if ok else INK, alpha=v)
        if k:
            t(ax, 14.9, yk, '✓ agree' if ok else 'Δ %.0e' % abs(P - lad[k - 1][1]), 11.5, TEAL if ok else MUT, alpha=v)
    t(ax, 11.0, 6.75 - 0.62 * len(lad) - 0.1, 'tolerance 1e-5;  converged value %.5f' % D['ladder_converged'], 12, MUT,
      alpha=seg(u, 0.94, 1.0))


# ------------------------------------------------------------------ diagram 4: strategy='auto'
AUTO_LEAVES = [(6.55, AMBER, 'constant density', 'one exact exponential per energy'),
               (5.05, TEAL, 'energy scan, one baseline', 'Magnus ladder, slabs shared by all energies'),
               (3.55, TEAL, 'baseline scan, one energy', 'cumulative scan: one pass along the path'),
               (2.05, TEAL, 'one point, modest phase', 'adaptive Magnus ladder'),
               (0.55, VIOLET, 'one point, extreme phase', 'hybrid: adiabatic transport + Magnus patches')]


def auto(fig, ax, u):
    t(ax, 0.8, 8.2, "STRATEGY='AUTO' PICKS THE SOLVER", 16, BLUE, DISP, weight=700, alpha=seg(u, 0, 0.08))
    v = seg(u, 0, 0.1)
    box(ax, 0.8, 3.4, 2.6, 1.2, ec=BLUE, lw=2, alpha=v)
    t(ax, 2.1, 4.0, 'your call', 15, INK, ha='center', va='center', alpha=v)
    for k, (y, col, cond, eng) in enumerate(AUTO_LEAVES):
        a = 0.1 + 0.16 * k
        v, w = seg(u, a + 0.04, a + 0.1), seg(u, a + 0.08, a + 0.14)
        arrow(ax, (3.45, 4.0), (5.2, y + 0.45), col=col, lw=1.4, p=seg(u, a, a + 0.06))
        box(ax, 5.25, y, 3.9, 0.9, ec=col, lw=1.4, alpha=v)
        t(ax, 5.45, y + 0.45, cond, 13, INK, va='center', alpha=v)
        arrow(ax, (9.2, y + 0.45), (9.75, y + 0.45), col=col, lw=1.4, alpha=w)
        t(ax, 9.85, y + 0.45, eng, 13, col, va='center', alpha=w)
    t(ax, 9.85, 0.15, 'certified, or handed back to the ladder', 11, MUT, alpha=seg(u, 0.9, 0.98))
    t(ax, 0.8, 2.4, 'strategy_info tells you\nwhich one answered', 12, MUT, linespacing=1.4, alpha=seg(u, 0.92, 1.0))


# ------------------------------------------------------------------ diagram 5: fast
def fast(fig, ax, u):
    D = load('diagrams.json')
    t(ax, 0.8, 8.2, 'THOUSANDS OF ENERGIES. ONE CALL.', 16, TEAL, DISP, weight=700, alpha=seg(u, 0, 0.08))
    E, P = np.asarray(D['fast']['E']), np.asarray(D['fast']['P'])
    s = seg(u, 0.1, 0.45)                                # the energies stream into the call
    for k, yy in enumerate(np.linspace(1.0, 7.4, 60)):
        ax.plot([0.9 + 1.8 * s, 2.2 + 1.26 * s], [yy, yy + (4.2 - yy) * s], color=BLUE, lw=1.1,
                alpha=(0.25 + 0.6 * (k % 7 == 0)) * (1 - seg(u, 0.38, 0.5)))
    t(ax, 1.55, 0.55, '2000 energies', 12, MUT, ha='center', alpha=seg(u, 0, 0.1) * (1 - seg(u, 0.4, 0.5)))
    v = seg(u, 0.02, 0.12)
    box(ax, 3.45, 3.4, 4.3, 1.6, fc='#0c1a14', ec=TEAL, lw=2.2 + 1.5 * np.exp(-((u - 0.47) / 0.04)**2), r=0.18, alpha=v)
    t(ax, 5.6, 4.45, 'osc_prob_3nu_earth(', 14, INK, ha='center', va='center', alpha=v)
    t(ax, 5.6, 3.95, '    energy, costhz=-0.8, ...)', 14, INK, ha='center', va='center', alpha=v)
    arrow(ax, (7.8, 4.2), (8.6, 4.2), col=TEAL, lw=2, alpha=seg(u, 0.46, 0.52))
    g = seg(u, 0.5, 0.95)
    a = fig.add_axes([0.56, 0.16, 0.40, 0.66])
    plain(a)
    k = int(g * len(E))
    if k > 1:
        a.plot(E[:k], P[:k], color=TEAL, lw=1.4)
    a.set_xscale('log')
    a.set_xlim(E[0], E[-1])
    a.minorticks_off()
    a.set_xticks([1, 3, 10, 30])
    a.set_xticklabels(['1', '3', '10', '30 GeV'] if g > 0 else [''] * 4)
    a.set_yticks([0, 0.2, 0.4])
    a.set_yticklabels(['0', '0.2', '0.4'] if g > 0 else [''] * 3)
    a.set_ylim(-0.01, max(0.45, P.max() * 1.05))
    monoticks(a)
    if g > 0:
        a.text(1.05, a.get_ylim()[1], r'$P(\nu_\mu \to \nu_e)$ through the Earth', color=INK, fontfamily=MONO,
               fontsize=13, va='top', alpha=g)


# ------------------------------------------------------------------ diagram 6: 2 to 5 flavors
def flavors(fig, ax, u):
    D = load('diagrams.json')
    t(ax, 0.8, 8.2, '2 TO 5 FLAVORS, READY-MADE', 16, AMBER, DISP, weight=700, alpha=seg(u, 0, 0.1))
    fcol = [BLUE, AMBER, TEAL, ROSE, VIOLET]
    names = [r'$\nu_e$', r'$\nu_\mu$', r'$\nu_\tau$', r'$\nu_{s1}$', r'$\nu_{s2}$']
    for g, dd in enumerate((2, 3, 4, 5)):
        U, gx, a = np.asarray(D['U'][str(dd)]), 0.9 + g * 3.85, 0.05 + 0.18 * g
        t(ax, gx, 7.25, {2: '2 flavors', 3: '3 flavors', 4: '3 + 1', 5: '3 + 2'}[dd], 15, INK, alpha=seg(u, a, a + 0.08))
        for k in range(dd):
            gr, yk, left = seg(u, a + 0.03 * k, a + 0.16 + 0.03 * k), 6.4 - k * 1.05, gx + 0.5
            for f in range(dd):
                w = 2.8 * U[f, k] * gr
                if w > 0:
                    ax.add_patch(Rectangle((left, yk), w, 0.6, fc=fcol[f], ec=BG, lw=0.8))
                left += w
            t(ax, gx + 0.4, yk + 0.3, r'$\nu_%d$' % (k + 1), 13, MUT, ha='right', va='center', alpha=gr)
    v = seg(u, 0.8, 0.95)
    for f in range(5):
        ax.add_patch(Rectangle((0.9 + f * 2.6, 0.55), 0.35, 0.35, fc=fcol[f], ec='none', alpha=v))
        t(ax, 1.35 + f * 2.6, 0.72, names[f] + ('' if f < 3 else ' sterile'), 13, MUT, va='center', alpha=v)
    t(ax, 0.9, 1.3, 'each bar: the flavor content of one mass state', 11.5, MUT, va='center', alpha=seg(u, 0.85, 1.0))


# name -> (draw function, beats, data it needs); beats match trailer.json
SCENES = {'opening': (opening, 14, 'opening'), 'switch': (switch, 12, 'adiabatic'),
          'burden': (burden, 10, None), 'hamiltonian': (hamiltonian, 8, 'diagrams'),
          'slabs': (slabs, 4, 'diagrams'), 'auto': (auto, 6, None), 'fast': (fast, 4, 'diagrams'),
          'flavors': (flavors, 3, 'diagrams')}
