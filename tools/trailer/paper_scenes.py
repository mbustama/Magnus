r"""The paper's figures, remade as trailer scenes: the trailer's fonts and colors, no paper
titles (the trailer's own title band names each scene), and elements uncovered progressively.
Each is drawn as a function of its progress u in [0, 1], like ``scenes.py``.  The numbers come
from ``build/paper.npz`` (``data.py paper``): the paper's own arrays where it cached them, the
rest computed with its settings.
"""
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, FancyArrowPatch, Polygon, Wedge

from common import BG, INK, MUT, BLUE, AMBER, TEAL, ROSE, VIOLET, MONO, DISP, seg
from scenes import load, plain, monoticks

# Probability, dark to bright: the trailer's ground, its blue, its teal, and a pale gold.
PCMAP = LinearSegmentedColormap.from_list('trailerP', ['#06080d', '#16305e', '#2f6fc0', '#3fd0a4', '#f3e7a1'])
DCMAP = LinearSegmentedColormap.from_list('trailerD', ['#5aa9ff', '#0f1626', '#f0a33c'])
R_E = 6371.0
PREM = [(0.0, 1221.5, '#f6d9a8', 'Inner core'), (1221.5, 3480.0, '#f0bd7e', 'Outer core'),
        (3480.0, 6346.6, '#c9895a', 'Mantle'), (6346.6, R_E, '#8f5a36', 'Crust')]


def P():
    return load('paper.npz')


def label(ax, x, y, s, size=13, col=MUT, **kw):
    ax.text(x, y, s, fontsize=size, color=col, fontfamily=MONO, **kw)


def dark(a):
    """A data panel on the trailer's ground: dark face, faint frame."""
    a.set_facecolor('#0b111c')
    for sp in a.spines.values():
        sp.set_visible(True)
        sp.set_color('#232c3f')
    return a


def axes(fig, rect, xlog=False):
    a = fig.add_axes(rect)
    plain(a)
    if xlog:
        a.set_xscale('log')
        a.minorticks_off()
    return a


_land = None


def land():
    """The continents of the paper's globes (notebook 28's LAND, [lon, lat] rings), read once."""
    global _land
    if _land is None:
        import json
        from common import BUILD
        _land = json.loads((BUILD / 'land.json').read_text())
    return _land


_land_images = {}


def land_image(detector, n=560):
    """The visible hemisphere as an image: ocean, land, and a mask for the coastlines.  Each pixel
    of the disc is taken back through the orthographic projection (centred 90 degrees south of the
    detector, as in the paper, so the detector sits at the top of the limb) to a latitude and
    longitude, and tested against notebook 28's land polygons.  Unlike clipping the polygons, this
    is right for any viewpoint, including continents that wrap around the limb."""
    if detector not in _land_images:
        from matplotlib.path import Path
        y, x = np.mgrid[1:-1:n * 1j, -1:1:n * 1j]
        rho = np.hypot(x, y)
        inside = rho <= 1.0
        c = np.arcsin(np.clip(rho, 0, 1))
        lat0, lon0 = np.radians(detector[0] - 90.0), np.radians(detector[1])
        with np.errstate(invalid='ignore', divide='ignore'):
            lat = np.arcsin(np.cos(c) * np.sin(lat0) + np.where(rho > 0, y * np.sin(c) * np.cos(lat0) / rho, 0))
            lon = lon0 + np.arctan2(x * np.sin(c), rho * np.cos(c) * np.cos(lat0) - y * np.sin(c) * np.sin(lat0))
        lon = (np.degrees(lon) + 180.0) % 360.0 - 180.0
        pts = np.column_stack([lon[inside], np.degrees(lat[inside])])
        is_land = np.zeros(len(pts), bool)
        for ring in land():
            path = Path(np.asarray(ring, float))
            for shift in (0.0, 360.0, -360.0):              # rings that cross the antimeridian
                is_land |= path.contains_points(pts + [shift, 0.0])
        mask = np.zeros(inside.shape)
        mask[inside] = is_land
        _land_images[detector] = (mask, inside)
    return _land_images[detector]


def earth_disc(a, cut=(270, 360), alpha=1.0, labels=True, detector=(42.5, 13.6)):
    """The Earth, seen with a quarter cut away to show the PREM layers, the continents on its face."""
    from matplotlib.colors import to_rgba
    mask, inside = land_image(detector)
    img = np.zeros(mask.shape + (4,))
    img[inside] = to_rgba('#15325a', alpha)
    img[mask > 0.5] = to_rgba('#3f6f4f', alpha)
    a.imshow(img, extent=(-1, 1, -1, 1), origin='upper', interpolation='bilinear', zorder=1)
    n = mask.shape[0]
    g = np.linspace(-1, 1, n)
    a.contour(g, g[::-1], mask, levels=[0.5], colors=['#6fa37e'], linewidths=0.7, alpha=alpha, zorder=2)
    for r0, r1, col, name in PREM[::-1]:
        a.add_patch(Wedge((0, 0), r1 / R_E, *cut, width=(r1 - r0) / R_E, fc=col, ec=BG, lw=0.8,
                          alpha=alpha, zorder=3))
    a.add_patch(Circle((0, 0), 1.0, fc='none', ec='#5aa9ff', lw=1.4, alpha=alpha, zorder=4))
    if labels:
        th = np.radians(0.5 * (cut[0] + cut[1]))
        for r0, r1, col, name in PREM[1:3]:             # outer core and mantle; the inner core is too small
            rm = 0.5 * (r0 + r1) / R_E
            a.text(rm * np.cos(th) + 0.02, rm * np.sin(th), name, fontsize=11, color=BG, fontfamily=MONO,
                   ha='center', va='center', rotation=np.degrees(th) + 90, alpha=alpha, zorder=5)
    a.set_xlim(-1.15, 1.15)
    a.set_ylim(-1.15, 1.15)
    a.set_aspect('equal')
    a.axis('off')


# ------------------------------------------------------------------ Earth oscillogram
def _oscillogram(a, cz, E, grid, k, cbar=None):
    img = np.full(grid.shape, np.nan)
    img[:k] = grid[:k]
    a.pcolormesh(cz, E, np.ma.masked_invalid(img).T, cmap=PCMAP, vmin=0, vmax=1, shading='gouraud',
                 rasterized=True)
    a.set_yscale('log')
    a.minorticks_off()
    a.set_xlim(cz[0], cz[-1])
    a.set_ylim(E[0], E[-1])
    a.set_facecolor('#0b111c')
    for s in a.spines.values():
        s.set_visible(True)
        s.set_color('#232c3f')


def earth_osc(fig, ax, u):
    d = P()
    cz, E, grid = d['osc_cz'], d['osc_e_gev'], d['osc_3nu']
    k = int(np.clip(seg(u, 0.05, 0.92), 0, 1) * len(cz))
    g = fig.add_axes([0.02, 0.10, 0.40, 0.72])
    earth_disc(g, cut=(270, 360), detector=(36.4, 137.3))             # Kamioka at the top
    c = cz[max(0, min(k, len(cz) - 1))]
    s = np.sqrt(max(0.0, 1 - c * c))
    end = (-2 * c * s, 1 - 2 * c * c)                  # detector at the top, looking down
    g.plot([0, end[0]], [1, end[1]], color=AMBER, lw=3, zorder=6, solid_capstyle='round')
    t = (u * 6) % 1                                      # a neutrino running up the chord
    g.scatter([end[0] * (1 - t)], [1 + (end[1] - 1) * (1 - t)], s=60, color='#fff', zorder=7)
    g.plot([0], [1.0], marker='v', ms=13, color=TEAL, zorder=8)
    g.text(0.06, 1.08, 'Detector', fontsize=13, color=TEAL, fontfamily=MONO, va='bottom')
    g.text(0, -1.13, r'$\cos\theta_z$ = %+.2f   chord %d km' % (c, d['osc_chord_km'][min(k, len(cz) - 1)]),
           fontsize=14, color=INK, fontfamily=MONO, ha='center', va='top')
    a = dark(fig.add_axes([0.49, 0.15, 0.40, 0.66]))
    _oscillogram(a, cz, E, grid, k)
    a.axvline(cz[min(k, len(cz) - 1)], color=AMBER, lw=2)
    a.axvline(-0.837, color='#ffffff', lw=0.8, ls=(0, (4, 3)), alpha=0.6)
    a.text(-0.83, E[-1] * 0.92, 'Core', color='#ffffff', fontsize=12, fontfamily=MONO, va='top', alpha=0.8)
    a.set_yticks([2, 5, 10, 20, 50])
    a.set_yticklabels(['2', '5', '10', '20', '50 GeV'])
    a.set_xticks([-1, -0.8, -0.6, -0.4, -0.2])
    a.tick_params(colors=MUT, labelsize=12, length=0)
    monoticks(a)
    a.set_xlabel(r'Direction, $\cos\theta_z$  (straight up through the Earth at $-1$)', color=MUT, fontsize=13)
    a.text(cz[0], E[-1] * 1.12, r'Survival probability, $P(\nu_\mu \to \nu_\mu)$, three flavors', color=INK,
           fontsize=15, fontfamily=MONO, va='bottom')
    cb = fig.add_axes([0.905, 0.15, 0.012, 0.66])
    cb.imshow(np.linspace(1, 0, 256)[:, None], aspect='auto', cmap=PCMAP, extent=(0, 1, 0, 1))
    cb.set_xticks([])
    cb.yaxis.tick_right()
    cb.set_yticks([0, 0.5, 1])
    cb.tick_params(colors=MUT, labelsize=12, length=0)
    monoticks(cb)
    for sp in cb.spines.values():
        sp.set_visible(False)


def earth_more(fig, ax, u):
    d = P()
    cz = d['osc_cz']
    panels = [('Three flavors', d['osc_e_gev'], d['osc_3nu'], 'GeV', [2, 10, 50], 0.0),
              ('3 + 1: one sterile neutrino', d['osc_e_tev'], d['osc_3p1'], 'TeV', [1, 3, 10, 30], 0.15),
              ('3 + 2: two sterile neutrinos', d['osc_e_tev'], d['osc_3p2'], 'TeV', [1, 3, 10, 30], 0.45)]
    for i, (name, E, grid, unit, ticks, a0) in enumerate(panels):
        v = 1.0 if i == 0 else seg(u, a0, a0 + 0.08)
        if v <= 0:
            continue
        k = len(cz) if i == 0 else int(seg(u, a0, a0 + 0.35) * len(cz))
        a = dark(fig.add_axes([0.05 + 0.315 * i, 0.16, 0.27, 0.62]))
        _oscillogram(a, cz, E, grid, k)
        a.set_yticks(ticks)
        a.set_yticklabels(['%g' % x for x in ticks[:-1]] + ['%g %s' % (ticks[-1], unit)])
        a.set_xticks([-1, -0.5])
        a.tick_params(colors=MUT, labelsize=12, length=0)
        monoticks(a)
        a.set_xlabel(r'$\cos\theta_z$', color=MUT, fontsize=13)
        a.text(cz[0], E[-1] * 1.1, name, color=[INK, ROSE, VIOLET][i], fontsize=16, fontfamily=MONO,
               va='bottom', alpha=v)


# ------------------------------------------------------------------ Fermilab to four sites
SITES = [('snolab', 'SNOLAB', BLUE), ('homestake', 'Homestake', TEAL), ('cern', 'CERN', AMBER),
         ('south_pole', 'South Pole', ROSE)]


def baselines(fig, ax, u):
    d = P()
    g = fig.add_axes([0.02, 0.12, 0.38, 0.70])
    earth_disc(g, cut=(200, 340), labels=False, detector=(41.8, -88.3))   # Fermilab at the top
    g.plot([0], [1], marker='o', ms=12, color='#ffffff', zorder=9)
    g.text(0, 1.07, 'Fermilab', fontsize=15, color=INK, fontfamily=MONO, ha='center', va='bottom')
    a = dark(fig.add_axes([0.47, 0.16, 0.48, 0.64]))
    a.set_xscale('log')
    a.minorticks_off()
    E = d['fnal_e_gev']
    for j, (key, name, col) in enumerate(SITES):
        a0 = 0.04 + 0.22 * j
        p = seg(u, a0, a0 + 0.2)
        if p <= 0:
            continue
        L = float(d['fnal_L_' + key])
        th = 2 * np.arcsin(min(1.0, L / (2 * R_E)))       # central angle of the chord
        side = [-1, 1, -1, 1][j]
        end = (side * np.sin(th), np.cos(th))
        g.plot([0, end[0] * p + 0 * (1 - p)], [1, 1 + (end[1] - 1) * p], color=col, lw=3, zorder=8,
               solid_capstyle='round')
        if p > 0.95:
            g.plot([end[0]], [end[1]], marker='o', ms=9, color=col, zorder=9)
        ax.text(0.6 + 3.3 * (j % 2), 0.95 - 0.55 * (j // 2), '%s, %s km' % (name, format(int(round(L)), ',')),
                fontsize=15, color=col, fontfamily=MONO, va='center', alpha=seg(u, a0, a0 + 0.05))
        k = int(p * len(E))
        if k > 1:
            a.plot(E[:k], d['fnal_' + key][:k], color=col, lw=2.6)
    a.set_xlim(E[0], E[-1])
    a.set_ylim(0, 0.4)
    a.set_xticks([0.3, 1, 3, 10])
    a.set_xticklabels(['0.3', '1', '3', '10 GeV'])
    a.set_yticks([0, 0.2, 0.4])
    a.tick_params(colors=MUT, labelsize=13, length=0)
    monoticks(a)
    a.text(E[0], 0.41, r'$P(\nu_\mu \to \nu_e)$ at the far detector', color=INK, fontsize=16, fontfamily=MONO,
           va='bottom')


# ------------------------------------------------------------------ CP violation
def cp(fig, ax, u):
    d = P()
    a = dark(fig.add_axes([0.20, 0.14, 0.40, 0.68]))
    dcp = d['cp_dcp']
    k = max(2, int(seg(u, 0.05, 0.85) * len(dcp)))
    for order, col, name in (('NO', BLUE, 'Normal ordering'), ('IO', AMBER, 'Inverted ordering')):
        x, y = d['cp_%s_nu' % order], d['cp_%s_nubar' % order]
        a.plot(x, y, color=col, lw=1, alpha=0.15)
        a.plot(x[:k], y[:k], color=col, lw=3)
        a.scatter([x[k - 1]], [y[k - 1]], s=90, color=col, zorder=5, edgecolor='#fff', lw=1.5)
        a.text(x[0], y[0], '', color=col)
    lo = min(d['cp_NO_nu'].min(), d['cp_NO_nubar'].min(), d['cp_IO_nu'].min(), d['cp_IO_nubar'].min())
    hi = max(d['cp_NO_nu'].max(), d['cp_NO_nubar'].max(), d['cp_IO_nu'].max(), d['cp_IO_nubar'].max())
    pad = 0.1 * (hi - lo)
    a.set_xlim(lo - pad, hi + pad)
    a.set_ylim(lo - pad, hi + pad)
    a.plot([0, 1], [0, 1], color=MUT, lw=0.8, ls=(0, (3, 3)))
    a.tick_params(colors=MUT, labelsize=12, length=0)
    monoticks(a)
    a.set_xlabel(r'$P(\nu_\mu \to \nu_e)$', color=MUT, fontsize=15)
    a.set_ylabel(r'$P(\bar\nu_\mu \to \bar\nu_e)$', color=MUT, fontsize=15)
    for s in ('left', 'bottom'):
        a.spines[s].set_visible(True)
        a.spines[s].set_color('#232c3f')
    ax.text(10.2, 6.4, r'$\delta_{CP}$ = %3d°' % round(np.degrees(dcp[k - 1])), fontsize=30, color=INK,
            fontfamily=MONO, va='center')
    ax.text(10.2, 5.3, 'Normal ordering', fontsize=16, color=BLUE, fontfamily=MONO, va='center')
    ax.text(10.2, 4.7, 'Inverted ordering', fontsize=16, color=AMBER, fontfamily=MONO, va='center')
    ax.text(10.2, 3.6, 'Fermilab → Homestake\n1,284 km, 2.5 GeV', fontsize=14, color=MUT, fontfamily=MONO,
            va='center', linespacing=1.5)


# ------------------------------------------------------------------ new physics
def bsm(fig, ax, u):
    d = P()
    a = dark(fig.add_axes([0.07, 0.16, 0.62, 0.64]))
    a.set_xscale('log')
    a.minorticks_off()
    E = d['bsm_e_gev']
    curves = [('bsm_std', INK, 'Standard', 0.02, (0, (5, 3)), 2.2),
              ('bsm_nsi', AMBER, '+ Non-standard interactions', 0.34, '-', 2.8),
              ('bsm_liv', VIOLET, '+ Lorentz violation', 0.64, '-', 2.8)]
    for key, col, name, a0, ls, lw in curves:
        p = seg(u, a0, a0 + 0.26)
        k = int(p * len(E))
        if k > 1:
            a.plot(E[:k], d[key][:k], color=col, lw=lw, ls=ls, zorder=4 if key == 'bsm_std' else 3)
        ax.text(11.5, {'bsm_std': 6.2, 'bsm_nsi': 5.4, 'bsm_liv': 4.6}[key], name, fontsize=15, color=col,
                fontfamily=MONO, va='center', alpha=seg(u, a0, a0 + 0.06))
    a.set_xlim(E[0], E[-1])
    a.set_ylim(0, 1.02)
    a.set_xticks([1, 2, 5, 10, 20, 40])
    a.set_xticklabels(['1', '2', '5', '10', '20', '40 GeV'])
    a.set_yticks([0, 0.5, 1])
    a.tick_params(colors=MUT, labelsize=13, length=0)
    monoticks(a)
    a.text(E[0], 1.05, r'$P(\nu_\mu \to \nu_\mu)$ through the Earth, 11,470 km', color=INK, fontsize=16,
           fontfamily=MONO, va='bottom')


# ------------------------------------------------------------------ flavor triangle
def _tri(f):
    f = np.atleast_2d(f)
    return f[:, 1] + 0.5 * f[:, 2], np.sqrt(3) / 2 * f[:, 2]


def triangle(fig, ax, u):
    d = P()
    a = dark(fig.add_axes([0.20, 0.08, 0.52, 0.78]))
    a.set_aspect('equal')
    a.axis('off')
    lo, hi = 0.30, 0.40
    corners = np.array([[hi, lo, lo], [lo, hi, lo], [lo, lo, hi]])
    cx, cy = _tri(corners)
    a.fill(cx, cy, fc='#0b111c', ec='#3a4660', lw=1.5)
    for v in np.arange(0.31, 0.40, 0.01):                  # gridlines of constant fraction
        for i in range(3):
            j, k = (i + 1) % 3, (i + 2) % 3
            p1 = np.zeros(3)
            p2 = np.zeros(3)
            p1[i] = p2[i] = v
            p1[j], p1[k] = lo, 1 - v - lo
            p2[j], p2[k] = 1 - v - lo, lo
            x, y = _tri(np.array([p1, p2]))
            a.plot(x, y, color='#1c2638', lw=0.8)
    for f, s, off in ((corners[0], r'more $\nu_e$', (0.0, -0.006)), (corners[1], r'more $\nu_\mu$', (0.003, -0.006)),
                      (corners[2], r'more $\nu_\tau$', (0.0, 0.004))):
        x, y = _tri(f)
        a.text(x[0] + off[0], y[0] + off[1], s, fontsize=15, color=MUT, fontfamily=MONO, ha='center',
               va='top' if off[1] < 0 else 'bottom')
    xv, yv = _tri(d['tri_vac'])
    a.scatter(xv, yv, s=140, color='#ffffff', zorder=6)
    a.text(xv[0] + 0.001, yv[0] - 0.003, 'Standard', fontsize=15, color=INK, fontfamily=MONO, ha='left', va='top')
    for key, col, name, a0 in (('tri_liv', AMBER, 'Lorentz violation', 0.12), ('tri_sterile', ROSE, 'Sterile neutrino', 0.52)):
        x, y = _tri(d[key])
        k = int(seg(u, a0, a0 + 0.36) * len(x))
        if k > 1:
            a.plot(x[:k], y[:k], color=col, lw=3.2, zorder=5, solid_capstyle='round')
            a.scatter([x[k - 1]], [y[k - 1]], s=70, color=col, zorder=7)
        if k == len(x):
            a.text(x[-1] + 0.002, y[-1], name, fontsize=15, color=col, fontfamily=MONO, va='center')
    a.set_xlim(cx.min() - 0.012, cx.max() + 0.016)
    a.set_ylim(cy.min() - 0.008, cy.max() + 0.005)
    ax.text(12.6, 6.2, 'Flavor at Earth,\nfrom a pion-decay\nsource 100 Mpc away', fontsize=14, color=MUT,
            fontfamily=MONO, va='center', linespacing=1.5)


# ------------------------------------------------------------------ the Sun, imaged
_sun_cache = {}


def _sun_image(i, n=560):
    if i not in _sun_cache:
        d = P()
        y, x = np.mgrid[-1:1:n * 1j, -1:1:n * 1j]
        r = np.hypot(x, y)
        img = np.interp(np.clip(r, 0, 1), np.linspace(0, 1, 2001), d['sun_P'][i])
        _sun_cache[i] = np.where(r <= 1, img, np.nan)
    return _sun_cache[i]


def sun(fig, ax, u):
    d = P()
    first = 2                                               # 30 GeV to 3 TeV: below and above, the Sun is uniform
    a = dark(fig.add_axes([0.24, 0.06, 0.44, 0.78]))
    a.set_aspect('equal')
    a.axis('off')
    scan = seg(u, 0.0, 0.25)                                  # the first image builds ray by ray
    E = d['sun_e_gev'][first:7]
    t = seg(u, 0.25, 0.97) * (len(E) - 1)
    i, f = int(np.floor(t)), t - np.floor(t)
    img = _sun_image(first + i)
    if f > 0 and i + 1 < len(E):
        img = (1 - f) * img + f * _sun_image(first + i + 1)
    if scan < 1:
        img = img.copy()
        img[:, int(scan * img.shape[1]):] = np.nan
    a.imshow(np.ma.masked_invalid(img), cmap=PCMAP, vmin=0, vmax=1, extent=(-1, 1, -1, 1), origin='lower',
             interpolation='bilinear')
    a.add_patch(Circle((0, 0), 1.0, fc='none', ec=AMBER, lw=1.8))
    if scan < 1:
        a.axvline(-1 + 2 * scan, color=AMBER, lw=2)
    a.set_xlim(-1.08, 1.08)
    a.set_ylim(-1.08, 1.08)
    e = E[min(len(E) - 1, int(round(t)))]
    txt = '%g MeV' % (e * 1e3) if e < 1 else ('%g GeV' % e if e < 1000 else '%g TeV' % (e / 1000))
    ax.text(12.3, 6.0, txt, fontsize=40, color=INK, fontfamily=[DISP, 'DejaVu Sans'], weight=700, va='center')
    ax.text(12.3, 5.0, r'$P(\nu_e \to \nu_e)$ of a diffuse', fontsize=15, color=MUT, fontfamily=MONO, va='center')
    ax.text(12.3, 4.5, 'flux crossing the Sun,', fontsize=15, color=MUT, fontfamily=MONO, va='center')
    ax.text(12.3, 4.0, 'one line of sight per pixel', fontsize=15, color=MUT, fontfamily=MONO, va='center')
    cb = fig.add_axes([0.07, 0.2, 0.012, 0.5])
    cb.imshow(np.linspace(1, 0, 256)[:, None], aspect='auto', cmap=PCMAP, extent=(0, 1, 0, 1))
    cb.set_xticks([])
    cb.set_yticks([0, 0.5, 1])
    cb.tick_params(colors=MUT, labelsize=12, length=0)
    monoticks(cb)
    for sp in cb.spines.values():
        sp.set_visible(False)


# ------------------------------------------------------------------ a buried body
def cavity(fig, ax, u):
    d = P()
    alpha, E, dP = d['cav_alpha'], d['cav_e_mev'], d['cav_dP']
    k = max(1, int(seg(u, 0.05, 0.9) * len(alpha)))
    g = fig.add_axes([0.04, 0.50, 0.92, 0.34])
    g.set_aspect('equal')
    g.axis('off')
    L0, R, D0 = 1500.0, 125.0, 750.0
    g.fill_between([-40, L0 + 60], -420, 420, color='#3a2618', alpha=0.55, lw=0)
    th = np.linspace(0, 2 * np.pi, 200)
    g.fill(D0 + R * np.cos(th), R * np.sin(th), fc=BLUE, alpha=0.35, ec=BLUE, lw=2)
    g.text(D0, -R - 25, 'Buried body, 10 g/cm³', fontsize=13, color=BLUE, fontfamily=MONO, ha='center', va='top')
    a_now = np.radians(alpha[k - 1])
    g.plot([0, L0 * np.cos(a_now)], [0, L0 * np.sin(a_now)], color=AMBER, lw=3, solid_capstyle='round')
    g.plot([0], [0], marker='o', ms=12, color='#ffffff')
    g.text(-20, -40, 'Source', fontsize=13, color=INK, fontfamily=MONO, ha='left', va='top')
    g.plot([L0 * np.cos(a_now)], [L0 * np.sin(a_now)], marker='s', ms=11, color=TEAL)
    g.text(L0 + 30, 0, 'Detector,\n1,500 km', fontsize=13, color=TEAL, fontfamily=MONO, va='center')
    g.set_xlim(-60, L0 + 260)
    g.set_ylim(-420, 420)
    a = dark(fig.add_axes([0.10, 0.06, 0.72, 0.38]))
    lim = np.abs(dP).max()
    img = np.full(dP.shape, np.nan)
    img[:, :k] = dP[:, :k]
    a.pcolormesh(E, alpha, np.ma.masked_invalid(img).T, cmap=DCMAP, vmin=-lim, vmax=lim, shading='gouraud',
                 rasterized=True)
    a.axhline(alpha[k - 1], color=AMBER, lw=1.5)
    a.set_xlim(E[0], E[-1])
    a.set_ylim(alpha[0], alpha[-1])
    a.set_facecolor('#0b111c')
    a.tick_params(colors=MUT, labelsize=12, length=0)
    monoticks(a)
    a.set_xlabel('Antineutrino energy [MeV]', color=MUT, fontsize=13)
    a.set_ylabel('Beam angle [°]', color=MUT, fontsize=13)
    ax.text(13.4, 2.6, 'Change in\n' + r'$P(\bar\nu_e \to \bar\nu_e)$' + '\ncaused by\nthe body', fontsize=14,
            color=MUT, fontfamily=MONO, va='center', linespacing=1.5)


# ------------------------------------------------------------------ geoneutrinos
GEO = [('local', 'Local crust, 100 km', TEAL, 10.0, 0.078), ('far_crust', 'Far crust, 3,400 km', BLUE, 20.0, -0.272),
       ('mantle', 'Mantle, 4,300 km', AMBER, 1000.0, -0.552), ('core', 'Through the core, 7,300 km', ROSE, 2800.0, -0.872)]


def geo(fig, ax, u):
    d = P()
    g = fig.add_axes([0.01, 0.10, 0.40, 0.74])
    earth_disc(g, cut=(0, 90), detector=(42.5, 13.6))                 # Gran Sasso at the top
    g.plot([0], [1.0], marker='*', ms=22, color='#ffffff', zorder=9)
    g.text(-0.05, 1.05, 'Detector', fontsize=14, color=INK, fontfamily=MONO, ha='right', va='bottom')
    Rdet = (R_E - 1.4) / R_E
    E = d['geo_e_mev']
    for j, (key, name, col, depth, c) in enumerate(GEO):
        a0 = 0.03 + 0.23 * j
        p = seg(u, a0, a0 + 0.18)
        if p <= 0:
            continue
        rs = (R_E - depth) / R_E
        L = float(d['geo_L_' + key]) / R_E                    # the true chord, from the source's depth
        if key == 'local':                                      # 100 km: drawn larger, or it is a dot
            xs, ys = 0.06, 0.995
        else:
            ca = (Rdet**2 + rs**2 - L**2) / (2 * Rdet * rs)
            ang = np.arccos(np.clip(ca, -1, 1))
            xs, ys = rs * np.sin(ang), rs * np.cos(ang)
        g.plot([xs, xs + (0 - xs) * p], [ys, ys + (Rdet - ys) * p], color=col, lw=3, zorder=8, solid_capstyle='round')
        g.plot([xs], [ys], marker='o', ms=10, color=col, zorder=9, mec=BG)
        a = dark(fig.add_axes([0.47, 0.64 - 0.165 * j, 0.50, 0.14]))
        plain(a)
        k = int(p * len(E))
        if k > 1:
            a.plot(E[:k], d['geo_' + key][:k], color=col, lw=0.6, rasterized=True)
        a.axhline(float(d['geo_vac_avg']), color=INK, lw=1.2, ls=(0, (4, 3)), alpha=0.7)
        a.set_xlim(1.8, 3.3)
        a.set_ylim(0, 1)
        a.set_yticks([0, 1])
        a.set_xticks([2, 2.5, 3] if j == 3 else [])
        if j == 3:
            a.set_xticklabels(['2', '2.5', '3 MeV'])
        a.tick_params(colors=MUT, labelsize=12, length=0)
        monoticks(a)
        a.text(1.82, 1.03, name, fontsize=13, color=col, fontfamily=MONO, va='bottom')
    ax.text(7.6, 8.2 - 0.0, r'$P(\bar\nu_e \to \bar\nu_e)$ from each production point', fontsize=14, color=INK,
            fontfamily=MONO, va='center', alpha=seg(u, 0.0, 0.1))


# ------------------------------------------------------------------ a jet inside a star
def jet(fig, ax, u):
    d = P()
    g = fig.add_axes([0.02, 0.18, 0.40, 0.62])
    g.set_aspect('equal')
    g.axis('off')
    for r, c in zip(np.linspace(1.0, 0.5, 7), ['#3a2a14', '#43301a', '#4c361e', '#553c22', '#5e4226', '#67482a', '#704e2e']):
        g.add_patch(Circle((0, 0), r, fc=c, ec='none'))
    g.add_patch(Circle((0, 0), 1.0, fc='none', ec=AMBER, lw=1.6))
    g.add_patch(Circle((0, 0), 0.5, fc='#8a6433', ec=AMBER, lw=1))
    g.add_patch(Polygon([[0, 0], [0.45, 0.14], [0.45, -0.14]], fc=AMBER, ec='none', alpha=0.9))
    g.text(0, 0.72, 'Hydrogen envelope', fontsize=12, color=INK, fontfamily=MONO, ha='center')
    g.text(-0.24, -0.1, 'Helium\ncore', fontsize=12, color=INK, fontfamily=MONO, ha='center', va='center')
    g.text(0.22, 0.2, 'Jet', fontsize=13, color=AMBER, fontfamily=MONO, ha='center')
    x = 0.46 + seg(u, 0.0, 0.5) * 1.1                        # the neutrino leaves through the envelope
    g.add_patch(FancyArrowPatch((0.46, 0), (1.62, 0), arrowstyle='-|>', mutation_scale=18, color=MUT, lw=1.2))
    g.scatter([x], [0], s=110, color='#ffffff', zorder=8)
    g.text(1.64, 0.1, 'To Earth', fontsize=13, color=MUT, fontfamily=MONO)
    g.set_xlim(-1.1, 2.1)
    g.set_ylim(-1.1, 1.1)
    a = dark(fig.add_axes([0.50, 0.18, 0.46, 0.62]))
    a.set_xscale('log')
    a.minorticks_off()
    E = d['jet_e_tev']
    for key, col, name, a0 in (('smooth', BLUE, 'Smooth envelope', 0.1), ('turbulent', TEAL, 'With turbulence', 0.38),
                               ('stepped', ROSE, 'With a sharp helium-core edge', 0.64)):
        k = int(seg(u, a0, a0 + 0.3) * len(E))
        if k > 1:
            a.plot(E[:k], d['jet_' + key][:k], color=col, lw=2.6)
        a.text(9e3, {'smooth': 0.22, 'turbulent': 0.16, 'stepped': 0.10}[key], name, fontsize=13, color=col, ha='right',
               fontfamily=MONO, alpha=seg(u, a0, a0 + 0.06))
    a.set_xlim(0.1, 1e4)
    a.set_ylim(0, 0.6)
    a.set_xticks([0.1, 1, 10, 100, 1000, 10000])
    a.set_xticklabels(['0.1', '1', '10', '100', '1000', '10⁴ TeV'])
    a.tick_params(colors=MUT, labelsize=12, length=0)
    monoticks(a)
    a.text(0.1, 0.62, r'$P(\nu_e \to \nu_e)$ at Earth', color=INK, fontsize=16, fontfamily=MONO, va='bottom')


# ------------------------------------------------------------------ a long-range force in the Sun
def lri(fig, ax, u):
    d = P()
    g = fig.add_axes([0.04, 0.50, 0.92, 0.34])
    g.set_aspect('equal')
    g.axis('off')
    s = seg(u, 0.0, 0.6)
    for cx, rng, col, name in ((0.0, None, INK, 'Standard'), (3.6, 0.1, BLUE, 'Range: a tenth of the Sun'),
                               (7.2, 1.0, ROSE, 'Range: the whole Sun')):
        g.add_patch(Circle((cx, 0), 1.0, fc='#3a2a14', ec=AMBER, lw=1.5))
        g.text(cx - 0.45, 0.55, 'Sun', fontsize=14, color=AMBER, fontfamily=MONO, ha='center', va='center')
        pos = cx + s * 1.25
        if rng is not None:
            g.add_patch(Circle((pos, 0), rng, fc=col if rng < 0.5 else 'none', ec=col, lw=2, alpha=0.5 if rng < 0.5 else 0.9))
        g.add_patch(FancyArrowPatch((cx, 0), (cx + 1.35, 0), arrowstyle='-|>', mutation_scale=14, color=MUT, lw=1))
        g.scatter([pos], [0], s=70, color='#ffffff', zorder=6)
        g.text(cx, -1.2, name, fontsize=13, color=col, fontfamily=MONO, ha='center', va='top')
    g.set_xlim(-1.4, 9.9)
    g.set_ylim(-1.6, 1.2)
    a = dark(fig.add_axes([0.10, 0.08, 0.80, 0.34]))
    a.set_xscale('log')
    a.minorticks_off()
    E = d['lr_e_mev']
    for key, col, ls in (('lr_std', INK, (0, (5, 3))), ('lr_01', BLUE, '-'), ('lr_1', ROSE, '-')):
        k = int(seg(u, 0.1, 0.9) * len(E))
        if k > 1:
            a.plot(E[:k], d[key][:k], color=col, lw=2.6, ls=ls)
    a.set_xlim(E[0], E[-1])
    lo = min(d['lr_std'].min(), d['lr_1'].min(), d['lr_01'].min())
    a.set_ylim(lo - 0.03, 0.58)
    a.set_xticks([0.1, 0.3, 1, 3, 10])
    a.set_xticklabels(['0.1', '0.3', '1', '3', '10 MeV'])
    a.tick_params(colors=MUT, labelsize=12, length=0)
    monoticks(a)
    a.text(E[0], 0.6, r'Average $P(\nu_e \to \nu_e)$ of solar neutrinos', color=INK, fontsize=15,
           fontfamily=MONO, va='bottom')


SCENES = {'earth_osc': (earth_osc, 10), 'earth_more': (earth_more, 6), 'baselines': (baselines, 8),
          'cp': (cp, 5), 'bsm': (bsm, 8), 'triangle': (triangle, 6), 'sun': (sun, 8), 'cavity': (cavity, 7),
          'geo': (geo, 8), 'jet': (jet, 7), 'lri': (lri, 7)}
