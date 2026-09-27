# The Magνs trailer

Generation scripts for the trailer planned in issue #99: a short advertisement for the
package, built from scenes computed with Magνs itself.  Every plot the trailer shows is
remade here, in the trailer's own look; no figure source in the repository is edited.

| File | What it is |
|---|---|
| `trailer.json` | The script and timeline: one entry per shot, in beats at 96 bpm, with its words, what moves, and where its picture comes from (a new `scene`, or the paper figure or notebook animation it will be remade from, `still_from`).  Also the decisions still open. |
| `common.py` | Paths, the look (colors, fonts), the two physics setups the new scenes use, and an independent reference integrator (no Magnus code) that checks the adiabatic scene. |
| `data.py` | Computes every number the new scenes show, with Magnus, and prints the checks the storyboard cites. |
| `scenes.py` | The new scenes, each drawn as a function of its progress `u` from 0 to 1: the opening's two tracks, the adiabatic switch, and six diagrams (the burden, the Hamiltonian, how Magnus works, how `strategy='auto'` picks a solver, many energies in one call, 2 to 5 flavors). |
| `render.py` | Renders each scene as a still (`u = 1`) and as frames and a clip. |
| `storyboard/` | `build.py` and its `template.html`: one page with a panel per shot, for review before the final render. |

## Running it

From the repository root, on one core each (`nice` keeps them out of the way of other work):

```bash
nice -n 19 python tools/trailer/data.py              # a few minutes; writes build/*.npz, build/diagrams.json
nice -n 19 python tools/trailer/render.py            # all scenes; or name some, or --stills-only
python tools/trailer/storyboard/build.py             # build/storyboard.html
```

Everything is written under `tools/trailer/build/`, which git ignores: the data, the stills
(1280×720), the PNG frames (30 fps) and the clips.  The PNG frames are the masters for the
final encode.  `render.py` writes H.264 MP4 when the `ffmpeg` on the PATH (or in `$FFMPEG`)
can; with the ffmpeg bundled with Playwright, which reads only JPEG and writes only VP8, it
writes a WebM for review instead.

The two fonts, IBM Plex Mono and Unbounded (both under the SIL Open Font License), are
downloaded from the Google Fonts repository on first use into `build/fonts/`.  The storyboard's
thumbnails of paper figures need PyMuPDF.

## What the scenes rest on

- **The opening**: P(νμ→νe) at 3 GeV over 0–10,000 km, in vacuum
  (`osc_prob_3nu_vacuum`) and through an illustrative density profile
  (`osc_prob_matter_std_potential`); NuFIT 6.1 parameters.
- **The code moment**: the user Hamiltonian shown on screen runs through
  `osc_prob_energy_baseline` without a warning and matches `osc_prob_matter_std_potential` to
  about 4e-8 (`data.py code`).
- **The adiabatic switch**: two flavors at 10 MeV crossing the MSW resonance slowly, through a
  shock, and slowly again.  One call per point under `strategy='auto'`: short paths are
  answered by the Magnus ladder and long ones by the hybrid, which finds the shock window
  itself; every checked point agrees with the independent reference to about 5e-5.  A sharper
  shock on a longer path is a case the hybrid gets wrong without a warning (issue #100); the
  scene stays out of that regime.
- **The diagrams**: the slab ladder is a real run at one point (fixed slab counts 2 to 21);
  the spectrum is a real 2000-energy `osc_prob_3nu_earth` call (no timing is claimed); the
  flavor bars are the package's own mixing matrices.
