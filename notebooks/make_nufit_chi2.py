# -*- coding: utf-8 -*-
r"""Extracts the one-dimensional Delta-chi^2 profiles from the NuFIT release data.

**The notebook data is generated.  Run this, do not hand-edit
``nufit_chi2.json``** -- and run it only when a new NuFIT release appears, since
it downloads some 60 MB from nu-fit.org that we deliberately do not keep.

Why this exists.  Notebook 26 asks how the *distribution* of oscillation
probabilities has moved as the global fits improved, which needs more than the
best-fit values :mod:`magnus.globaldefs` already ships -- it needs the
likelihood.  NuFIT publishes one per release as a set of Delta-chi^2
projections.

What is downloaded, and what is kept.  Each release's normal-ordering file is
4--5 MB compressed and about 1.5 million lines, almost all of which is a single
three-dimensional projection this notebook does not use.  Only the **six
one-dimensional projections** are extracted -- they sit in the last ~900 lines
of every file, are present in every release from v2.0 onward, and come to about
20 kB per release.  That is small enough to commit, which is the point: the
notebook must rebuild from the repository alone, without the network.

**Attribution.**  The data are the NuFIT collaboration's, redistributed here in
extract.  Cite the corresponding NuFIT paper and http://www.nu-fit.org/ when
using them.  Each entry records the URL it came from.

Two limits worth knowing, both stated in notebook 26 as well:

1. **v1.0 through v1.3 publish no machine-readable chi^2** -- only figures
   (``vNN.fig-chisq-glob.jpg``).  They are absent here rather than digitised,
   because reading numbers off a raster plot produces values that look
   authoritative and are not.
2. **These are one-dimensional marginals.**  Sampling each parameter
   independently from its own profile ignores the correlations, and the
   dCP--theta23 correlation in particular is real.  The pairwise and (from v5.0)
   three-dimensional projections are in the same source files if that matters
   for what you are doing; ``SECTIONS_2D`` names the one to reach for first.
"""

import gzip
import json
import lzma
import pathlib
import re
import urllib.request
import zlib


HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE/'nufit_chi2.json'
BASE = 'http://www.nu-fit.org/sites/default/files/'

# Normal ordering, and the with-SK-atmospheric variant wherever the release
# splits on it.  v6.0 renamed that category from 'SKyes' to 'TByes'; v2.1 offers
# LEM/LID instead and we take LEM, its primary.  v1.0-v1.3 are absent: figures
# only, see the module docstring.
RELEASES = {
    'NuFIT 2.0': 'v20.release-data-NO.txt.gz',
    'NuFIT 2.1': 'v21.release-LEM-NO.txt.gz',
    'NuFIT 2.2': 'v22.release-data-NO.txt.gz',
    'NuFIT 3.0': 'v30.release-data-NO.txt.gz',
    'NuFIT 3.1': 'v31.release-data-NO.txt.gz',
    'NuFIT 3.2': 'v32.release-data-NO.txt.gz',
    'NuFIT 4.0': 'v40.release-SKyes-NO.txt.gz',
    'NuFIT 4.1': 'v41.release-SKyes-NO.txt.gz',
    'NuFIT 5.0': 'v50.release-SKyes-NO.txt.xz',
    'NuFIT 5.1': 'v51.release-SKyes-NO.txt.xz',
    'NuFIT 5.2': 'v52.release-SKyes-NO.txt.xz',
    'NuFIT 5.3': 'v53.release-SKyes-NO.txt.xz',
    'NuFIT 6.0': 'v60.release-TByes-NO.txt.xz',
    'NuFIT 6.1': 'v61.release-TByes-NO.txt.xz',
}

# The six one-dimensional projections, keyed by the name used in the file and
# mapped to the Magnus parameter the column describes.  The values are the
# *file's* variable, not the Magnus one: T12 is sin^2(theta12) and DMS is
# log10(Delta_m21^2), and notebook 26 converts.
SECTIONS_1D = {
    'T12': 'sin^2(theta12)',
    'T13': 'sin^2(theta13)',
    'T23': 'sin^2(theta23)',
    'DCP': 'Delta_CP/deg',
    'DMS': 'log10(Delta_m21^2/eV^2)',
    'DMA': 'Delta_m31^2/(1e-3 eV^2)',
}

# The pairwise projection to reach for first if the independence assumption
# above is not good enough for what you are doing.  Not extracted by default.
SECTIONS_2D = 'T23/DCP'


# ----------------------------------------------------------------------
# v1.0 - v1.3: Gaussian approximations from the published parameter tables
# ----------------------------------------------------------------------
#
# These four releases publish no machine-readable chi^2, so the profiles above
# do not exist for them.  What they do publish is a parameter table --
# ``vNN.tbl-parameters.pdf`` -- and the entries below are transcribed from the
# **Free Fluxes + RSBL**, normal-ordering column of each.
#
# **Widths come from the 3-sigma range, not from the quoted +-1sigma.**  That is
# deliberate, and it is the whole reason this table stores ranges rather than
# sigmas.  NuFIT brackets a *local* minimum, and the +-1sigma printed inside the
# brackets is the local curvature: v1.3 quotes Delta_m31^2 = [2.458 +- 0.002]
# against a 3-sigma range of +-0.137, so believing the printed 1sigma would make
# that parameter some twenty times too sharp.  Deriving sigma from the 3-sigma
# range instead reproduces the published interval by construction and is immune
# to the bracket convention.
#
# The asymmetry is kept: sigma_lo = (bfp - lo)/3 and sigma_hi = (hi - bfp)/3
# define a two-piece ("split") normal, which the notebook samples directly.
#
# Where a release quotes two disconnected sin^2(theta23) minima joined by
# ``(+)``, the **unbracketed** one is the global best fit, and that is the
# centre used here -- verified against the values
# :func:`magnus.globaldefs.load_nufit_params` ships for the same release.
#
# delta_CP is unconstrained at 3 sigma in all four (the range is the full
# circle), so it is marked 'uniform' rather than given a spurious width.
#
# Transcribed 2026-08-09 from the PDFs at
# http://www.nu-fit.org/sites/default/files/vNN.tbl-parameters.pdf
GAUSSIAN_TABLES = {
    #                     bfp,     3sigma_lo, 3sigma_hi
    'NuFIT 1.0': {'url': BASE + 'v10.tbl-parameters.pdf', 'year': 2012, 'params': {
        'T12': (0.302, 0.267, 0.344),
        'T23': (0.413, 0.342, 0.667),
        'T13': (0.0227, 0.0156, 0.0299),
        'DCP': (300.0, 0.0, 360.0),
        'DMS': (7.50, 7.00, 8.09),
        'DMA': (2.473, 2.276, 2.695)}},
    'NuFIT 1.1': {'url': BASE + 'v11.tbl-parameters.pdf', 'year': 2013, 'params': {
        'T12': (0.306, 0.271, 0.346),
        'T23': (0.437, 0.357, 0.654),
        'T13': (0.0231, 0.0161, 0.0299),
        'DCP': (341.0, 0.0, 360.0),
        'DMS': (7.45, 6.98, 8.05),
        'DMA': (2.421, 2.248, 2.612)}},
    'NuFIT 1.2': {'url': BASE + 'v12.tbl-parameters.pdf', 'year': 2013, 'params': {
        'T12': (0.306, 0.271, 0.346),
        'T23': (0.593, 0.366, 0.663),
        'T13': (0.0231, 0.0173, 0.0288),
        'DCP': (266.0, 0.0, 360.0),
        'DMS': (7.45, 6.98, 8.05),
        'DMA': (2.417, 2.247, 2.623)}},
    'NuFIT 1.3': {'url': BASE + 'v13.tbl-parameters.pdf', 'year': 2014, 'params': {
        'T12': (0.304, 0.270, 0.344),
        'T23': (0.577, 0.385, 0.644),
        'T13': (0.0219, 0.0188, 0.0251),
        'DCP': (251.0, 0.0, 360.0),
        'DMS': (7.50, 7.03, 8.09),
        'DMA': (2.458, 2.325, 2.599)}},
}

# Units in which the table above is written, matched to the profile files so the
# notebook can treat both kinds of release identically.  DMS is the one that
# differs: the chi^2 files carry log10(Delta_m21^2/eV^2), the tables carry
# Delta_m21^2 in 1e-5 eV^2.
GAUSSIAN_UNITS = {'T12': 'sin^2(theta12)', 'T13': 'sin^2(theta13)',
                  'T23': 'sin^2(theta23)', 'DCP': 'Delta_CP/deg',
                  'DMS': 'Delta_m21^2/(1e-5 eV^2)',
                  'DMA': 'Delta_m31^2/(1e-3 eV^2)'}


# The second sin^2(theta23) minimum each release quotes alongside its best fit,
# joined by "(+)" in the tables.  v1.1 quotes only one, so it is absent here.
# Used to check that the extracted curve reproduces *both* minima -- the whole
# reason these releases are worth reading off the figures rather than
# approximating by a Gaussian, which has only one.
T23_SECOND_MINIMUM = {
    'NuFIT 1.0': 0.594,
    'NuFIT 1.2': 0.446,
    'NuFIT 1.3': 0.451,
}


def _pdf_content(blob: bytes) -> str:
    r"""The largest decompressed stream in a PDF -- its page content."""
    best = b''
    for chunk in re.findall(rb'stream\r?\n(.*?)endstream', blob, re.S):
        try:
            dec = zlib.decompress(chunk)
        except Exception:
            continue
        if len(dec) > len(best):
            best = dec
    return best.decode('latin-1')


def _polylines(txt: str) -> list:
    r"""Every stroked polyline in a page, as ``(rgb, solid, points)``.

    This is what makes the v1.x figures usable at all: they are **vector** PDFs,
    so the curves are stroke operators rather than pixels.  The colour operator
    separates the orderings exactly (``1 0 0 RG`` is NO, ``0 0 1 RG`` is IO) and
    the dash-array operator separates the analyses exactly (solid is
    Free+RSBL, dashed is Huber).  No digitising, and no guessing which curve is
    which where they overlap.
    """
    rgb, solid, cur, out = (0, 0, 0), True, [], []
    tokens = re.findall(r'\[[^\]]*\]\s*\d*\.?\d*\s*d'
                        r'|[-0-9.]+\s+[-0-9.]+\s+[-0-9.]+\s+RG'
                        r'|[-0-9.]+\s+[-0-9.]+\s+[ml]|[Sf]\b', txt)
    for tok in tokens:
        if tok.endswith('RG'):
            rgb = tuple(round(float(v)) for v in tok.split()[:3])
        elif tok.endswith('d'):
            solid = tok[tok.index('[') + 1:tok.index(']')].strip() == ''
        elif tok.endswith('m'):
            if len(cur) > 1:
                out.append((rgb, solid, cur))
            cur = [[float(tok.split()[0]), float(tok.split()[1])]]
        elif tok.endswith('l'):
            cur.append([float(tok.split()[0]), float(tok.split()[1])])
        elif tok in ('S', 'f'):
            if len(cur) > 1:
                out.append((rgb, solid, cur))
            cur = []
    if len(cur) > 1:
        out.append((rgb, solid, cur))
    return out


def _panel_grid(lines: list) -> tuple:
    r"""The 3x2 arrangement of plot boxes, from the long black frame edges."""
    black = [p for rgb, s, p in lines if rgb == (0, 0, 0) and s]
    def span(p, i):
        vals = [q[i] for q in p]
        return max(vals) - min(vals)
    def mean(p, i):
        return round(sum(q[i] for q in p)/len(p))
    horiz = [p for p in black if span(p, 1) < 5 and span(p, 0) > 800]
    vert = [p for p in black if span(p, 0) < 5 and span(p, 1) > 800]
    ys = sorted({mean(p, 1) for p in horiz})
    xs = sorted({mean(p, 0) for p in vert})
    return ([(ys[i], ys[i + 1]) for i in range(0, len(ys), 2)],
            [(xs[i], xs[i + 1]) for i in range(0, len(xs), 2)])


def _red_curve(lines: list, box: tuple):
    r"""The longest solid-red (NO, Free+RSBL) polyline inside a plot box."""
    xl, xr, yb, yt = box
    inside = []
    for rgb, s, p in lines:
        if rgb != (1, 0, 0) or not s or len(p) <= 3:
            continue
        xsv = [q[0] for q in p]
        ysv = [q[1] for q in p]
        if min(xsv) >= xl - 5 and max(xsv) <= xr + 5 \
                and min(ysv) >= yb - 5 and max(ysv) <= yt + 5:
            inside.append(p)
    return max(inside, key=len) if inside else None


def _figure_panels(blob: bytes) -> dict:
    r"""Maps each plot box to its parameter, detecting which layout this release uses.

    The arrangement is not the same in every release -- v1.0 and v1.1 put
    Delta_m21^2 top-left, v1.2 and v1.3 put sin^2(theta12) there.  Rather than
    hard-code a per-release table, the layout is *detected*: the Delta_m^2 panel
    has a split axis (negative for IO, positive for NO), so its NO curve fills
    only about a quarter of its box while every other panel's fills most of it.
    """
    lines = _polylines(_pdf_content(blob))
    rows, cols = _panel_grid(lines)
    yb, yt = rows[1]
    fill = []
    for xl, xr in cols:
        p = _red_curve(lines, (xl, xr, yb, yt))
        if p is None:
            fill.append(9.9)
        else:
            xsv = [q[0] for q in p]
            fill.append((max(xsv) - min(xsv))/(xr - xl))
    dma_left = fill[0] < fill[1]
    layout = ({(2, 0): 'DMS', (2, 1): 'T12', (1, 0): 'DMA', (1, 1): 'T23'} if dma_left
              else {(2, 0): 'T12', (2, 1): 'DMS', (1, 0): 'T23', (1, 1): 'DMA'})
    layout.update({(0, 0): 'T13', (0, 1): 'DCP'})

    out = {}
    for r, (yb, yt) in enumerate(rows):
        for c, (xl, xr) in enumerate(cols):
            p = _red_curve(lines, (xl, xr, yb, yt))
            if p is not None:
                out[layout[(r, c)]] = {'box': (xl, xr, yb, yt), 'pts': p}
    return out


CHI2_AXIS_TOP = 15.0
r"""float: the Delta_chi^2 shown at the top of every panel in the v1.x figures."""


def _to_profile(entry: dict, lo: float, hi: float, is_dcp: bool) -> dict:
    r"""One panel -> a (x, Delta_chi^2) profile in physical units.

    The vertical axis is read straight off the frame, which runs 0 to 15 in
    every panel.  The horizontal axis is anchored on the **published 3-sigma
    range**: the two outermost crossings of Delta_chi^2 = 9 are mapped onto it.
    That needs no axis-label parsing and reproduces the published interval by
    construction, leaving the *minimum* free as an independent check.

    delta_CP is the exception -- it is unconstrained at 3 sigma in all four of
    these releases, so there are no crossings to anchor on and the frame itself
    is mapped to [0, 360).
    """
    xl, xr, yb, yt = entry['box']
    pts = sorted(entry['pts'], key=lambda q: q[0])
    xp = [q[0] for q in pts]
    chi2 = [CHI2_AXIS_TOP*(q[1] - yb)/(yt - yb) for q in pts]

    if is_dcp:
        x = [lo + (hi - lo)*(v - xl)/(xr - xl) for v in xp]
    else:
        hits = []
        for i in range(len(xp) - 1):
            a, b = chi2[i], chi2[i + 1]
            if (a - 9.0)*(b - 9.0) < 0:
                hits.append(xp[i] + (9.0 - a)/(b - a)*(xp[i + 1] - xp[i]))
        if len(hits) < 2:
            raise SystemExit('no Delta_chi^2 = 9 crossings to anchor on')
        scale = (hi - lo)/(max(hits) - min(hits))
        x = [lo + (v - min(hits))*scale for v in xp]
    return {'x': x, 'chi2': chi2}


def _tolerance(x: list, lo: float, hi: float) -> float:
    r"""How precisely a minimum can be located on this particular polyline.

    The curves are stored at whatever resolution the figure was drawn with --
    as few as fifteen vertices across a panel -- so the position of a minimum
    is only meaningful to about the vertex spacing.  Demanding better than that
    fails a correct extraction for being coarsely sampled, which is what it did
    on the v1.1 sin^2(theta13) panel: right to within one vertex, rejected by a
    flat five-percent rule.
    """
    if len(x) < 2:
        return 0.05*(hi - lo)
    gaps = sorted(abs(x[i + 1] - x[i]) for i in range(len(x) - 1))
    spacing = gaps[len(gaps)//2]
    return max(0.05*(hi - lo), 1.5*spacing)


def build_figures() -> dict:
    r"""v1.0 - v1.3, read out of the vector figures and validated.

    Two checks, both against numbers this function never uses as input:

    * for the five constrained parameters, the curve's minimum must land on the
      published best fit;
    * for delta_CP, the curve must never reach Delta_chi^2 = 9, which is what
      the published "0 -> 360" 3-sigma range asserts.

    sin^2(theta23) is exempt from the first check wherever the release quotes
    two minima: the curve legitimately has two, and demanding a single one was
    the wrong question -- the pair is verified against the published pair
    instead.
    """
    out = {}
    for release, block in GAUSSIAN_TABLES.items():
        url = block['url'].replace('.tbl-parameters.pdf', '.fig-chisq-glob.pdf')
        print('  %-12s %s' % (release, url.rsplit('/', 1)[-1]), flush=True)
        with urllib.request.urlopen(url, timeout=600) as response:
            panels = _figure_panels(response.read())

        profiles, checks = {}, []
        for name, (bfp, lo, hi) in block['params'].items():
            is_dcp = (name == 'DCP')
            profiles[name] = _to_profile(panels[name], lo, hi, is_dcp)
            x, chi2 = profiles[name]['x'], profiles[name]['chi2']
            if is_dcp:
                ok = max(chi2) < 9.0
            elif name == 'T23' and release in T23_SECOND_MINIMUM:
                second = T23_SECOND_MINIMUM[release]
                tol = _tolerance(x, lo, hi)
                minima = [x[i] for i in range(1, len(x) - 1)
                          if chi2[i] <= chi2[i - 1] and chi2[i] <= chi2[i + 1]]
                ok = any(abs(m - bfp) < tol for m in minima) and \
                     any(abs(m - second) < tol for m in minima)
            else:
                ok = abs(x[chi2.index(min(chi2))] - bfp) < _tolerance(x, lo, hi)
            checks.append((name, ok))

        failed = [n for n, ok in checks if not ok]
        if failed:
            raise SystemExit('%s: figure extraction failed validation for %s'
                             % (release, ', '.join(failed)))
        out[release] = {'url': url, 'year': block['year'], 'profiles': profiles,
                        'units': dict(GAUSSIAN_UNITS)}
    return out


def build_gaussians() -> dict:
    r"""Turns each published 3-sigma range into a two-piece normal.

    ``delta_CP`` is flagged uniform wherever its range covers the full circle,
    which in these four releases it always does.
    """
    out = {}
    for release, block in GAUSSIAN_TABLES.items():
        entry = {'url': block['url'], 'year': block['year'], 'params': {}}
        for name, (bfp, lo, hi) in block['params'].items():
            if name == 'DCP' and lo <= 0.0 and hi >= 360.0:
                entry['params'][name] = {'kind': 'uniform', 'lo': 0.0, 'hi': 360.0,
                                         'bfp': bfp}
                continue
            if not lo < bfp < hi:
                raise SystemExit('%s/%s: best fit %g outside its 3-sigma range '
                                 '[%g, %g]' % (release, name, bfp, lo, hi))
            entry['params'][name] = {'kind': 'split_normal', 'bfp': bfp,
                                     'sigma_lo': (bfp - lo)/3.0,
                                     'sigma_hi': (hi - bfp)/3.0,
                                     'lo': lo, 'hi': hi}
        out[release] = entry
    return out


def read_release(url: str) -> list:
    r"""Downloads one release file and returns its lines.

    Held in memory rather than written to disk: the extract is the artefact
    worth keeping, and a stale multi-megabyte download in the working tree is
    the kind of thing that quietly gets committed.
    """
    with urllib.request.urlopen(url, timeout=600) as response:
        blob = response.read()
    opener = lzma.decompress if url.endswith('.xz') else gzip.decompress
    return opener(blob).decode('utf-8', errors='replace').splitlines()


def extract_1d(lines: list) -> dict:
    r"""Pulls the six one-dimensional profiles out of a release file.

    A section runs from its ``# NAME projection:`` header to the next header or
    end of file.  The two-dimensional sections share the ``A/B projection``
    spelling, so the split on ``/`` is what distinguishes them -- matching on
    the name alone would collect ``T23/DCP`` under ``T23``.
    """
    out, current = {}, None
    for line in lines:
        stripped = line.strip()
        if stripped.startswith('#'):
            current = None
            if 'projection:' in stripped:
                name = stripped.lstrip('#').split('projection:')[0].strip()
                if name in SECTIONS_1D:          # excludes 'T23/DCP' and friends
                    current = name
                    out[name] = {'x': [], 'chi2': []}
            continue
        if current is None or not stripped:
            continue
        parts = stripped.split()
        if len(parts) != 2:
            # A one-dimensional section must have exactly two columns.  Anything
            # else means the section boundaries moved, and guessing would put
            # the wrong numbers under the right key.
            raise ValueError('%s projection: expected 2 columns, got %d in %r'
                             % (current, len(parts), stripped))
        out[current]['x'].append(float(parts[0]))
        out[current]['chi2'].append(float(parts[1]))
    return out


def build():
    data = {'_source': 'http://www.nu-fit.org/',
            '_attribution': ('Delta-chi^2 profiles from the NuFIT collaboration; '
                             'cite the corresponding NuFIT paper and nu-fit.org. '
                             'Extracted by notebooks/make_nufit_chi2.py.'),
            '_ordering': 'NO',
            '_note': ('One-dimensional marginals only.  Sampling them '
                      'independently ignores the dCP-theta23 correlation; see '
                      'the module docstring.'),
            'releases': {}}

    print('  reading the v1.x figures (vector PDFs)...', flush=True)
    for release, block in build_figures().items():
        data['releases'][release] = {'url': block['url'],
                                     'profiles': block['profiles'],
                                     'source': 'figure'}

    for release, filename in RELEASES.items():
        url = BASE + filename
        print('  %-12s %s' % (release, filename), flush=True)
        profiles = extract_1d(read_release(url))
        missing = set(SECTIONS_1D) - set(profiles)
        if missing:
            raise SystemExit('%s: missing projections %s' % (release, sorted(missing)))
        for name, profile in profiles.items():
            if len(profile['x']) < 10:
                raise SystemExit('%s/%s: only %d points'
                                 % (release, name, len(profile['x'])))
        data['releases'][release] = {'url': url, 'profiles': profiles,
                                     'source': 'chi2-file'}

    OUT.write_text(json.dumps(data, separators=(',', ':')))
    total = sum(len(p['x']) for r in data['releases'].values()
                for p in r['profiles'].values())
    print('  wrote %s (%.0f kB, %d releases, %d points)'
          % (OUT.name, OUT.stat().st_size/1024, len(data['releases']), total))


if __name__ == '__main__':
    build()
