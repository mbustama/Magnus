# The Magνs trailer

Generation scripts for the trailer planned in issue #99: a short advertisement for the
package, built from scenes computed with Magνs itself.  Every plot the trailer shows is
remade here, in the trailer's own look; no figure source in the repository is edited.

| File | What it is |
|---|---|
| `trailer.json` | The script and timeline: one entry per shot, in beats at 96 bpm, with its words, what moves, and where its picture comes from (a new `scene`, or the paper figure or notebook animation it will be remade from, `still_from`).  Also the decisions still open. |
| `common.py` | Paths, the look (colors, fonts), the two physics setups the new scenes use, and an independent reference integrator (no Magnus code) that checks the adiabatic scene. |
| `data.py` | Computes every number the new scenes show, with Magnus, and prints the checks the storyboard cites. |
| `scenes.py` | The new scenes, each drawn as a function of its progress `u` from 0 to 1. |
| `render.py` | Renders each scene as a still (`u = 1`), as PNG frames, and as a clip. |
| `storyboard/` | `build.py` and its `template.html`: one page with a panel per shot, for review. |
| `cut.py` | Assembles the whole trailer, shot by shot, as 1920×1080 frames, and muxes them with the music into `build/magnus_trailer.mp4`. |
| `music.py` | The music: an original piece synthesized from sine waves and noise (nothing sampled, nothing to license), following each act's cue; the climax melody plays the Earth spectrum of the "Fast" diagram. |

## Creating the animations

### 1. What you need

- The repository's own environment (`pip install -e .` from the root): NumPy, SciPy and
  Matplotlib (the scenes were made with 3.10.8), which brings Pillow.
- `ffmpeg` on the PATH, or its path in `$FFMPEG`, to turn frames into clips.
- For the storyboard only: PyMuPDF (`pip install pymupdf`), for thumbnails of paper figures.
- Internet access on the first run: the two fonts, IBM Plex Mono and Unbounded (both under the
  SIL Open Font License), are downloaded from the Google Fonts repository into
  `build/fonts/`.  After that, nothing is fetched.

Run everything from the repository root.  `nice -n 19` keeps a render out of the way of other
work; each script uses one core, so scenes can be rendered in parallel from separate shells.

### 2. Compute the data

```bash
nice -n 19 python tools/trailer/data.py
```

This writes `build/opening.npz`, `build/adiabatic.npz` and `build/diagrams.json`, and prints
one check line per step.  They should read:

```
opening: max P vacuum 0.108, matter 0.233
code moment: max |P - wrapper| = 3.7e-08, warnings: none
adiabatic: ... engines {'magnus': 15, 'hybrid': 225}, warnings none, max |P - reference| at 10 points = 5.2e-05, windows [(44993.6, 45005.1)]
diagrams: ladder [(2, 0.01118), ..., (21, 0.10203)], converged 0.10202; ...; |U|^2 rows sum to 1: True
```

A step can be run on its own: `data.py opening`, `code`, `adiabatic` or `diagrams`.  The
adiabatic step is the slow one (about a minute here): it makes 240 separate calls, so that
each point is answered by the engine `strategy='auto'` picks for it.

### 3. Render the scenes

```bash
nice -n 19 python tools/trailer/render.py                       # every scene, 1280x720 (review)
nice -n 19 python tools/trailer/render.py switch opening        # only these scenes
nice -n 19 python tools/trailer/render.py --stills-only         # the finished frames only
nice -n 19 python tools/trailer/render.py --width 1920          # final size, 1920x1080
```

The scenes, with their length at 96 bpm (the same as their shot in `trailer.json`):

| Scene | Beats | Seconds | Frames | What moves |
|---|---|---|---|---|
| `opening` | 14 | 8.8 | 262 | One neutrino rides two tracks, vacuum and matter; each P(νμ→νe) draws itself up to it |
| `burden` | 10 | 6.3 | 188 | 25 solver files fill an experiments × theories grid, then converge on one H(E, x) |
| `hamiltonian` | 8 | 5.0 | 150 | H = H_vac/E + V(x)·P_e + H_new(x) writes itself; each term lights its card |
| `auto` | 6 | 3.8 | 112 | The five branches of `strategy='auto'` light in turn |
| `switch` | 12 | 7.5 | 225 | The neutrino crosses three resonances; the shock window turns amber and is magnified; P_ee lands one dot per call, with the engine that answered |
| `flavors` | 3 | 1.9 | 56 | Flavor-content bars grow from 2 to 5 flavors |
| `slabs` | 4 | 2.5 | 75 | The path is sliced, the exponentials multiply, the real ladder converges |
| `fast` | 4 | 2.5 | 75 | 2000 energies stream into one call; the spectrum draws |

For each scene the output is:

- `build/stills/<scene>.png`: the finished frame;
- `build/frames/<scene>/f0000.png ...`: every frame, the masters for the final edit (the last
  12 % of each scene holds the finished frame);
- `build/clips/<scene>.mp4`: H.264, when the ffmpeg can read PNG frames and has `libx264`.
  Otherwise (the ffmpeg bundled with Playwright reads only JPEG and writes only VP8) the
  frames go through quality-95 JPEG into `build/clips/<scene>.webm`, good enough for review.

To change a scene, edit its function in `scenes.py` (it receives the figure, a full-frame
axes in a 16 × 9 coordinate box, and `u`), then re-render just that scene.  `seg(u, a, b)`
from `common.py` gives the eased progress of a step that runs from `u = a` to `u = b`.

### 4. Review them in the storyboard

```bash
python tools/trailer/storyboard/build.py        # -> build/storyboard.html
```

One self-contained page: the timeline, then a panel per shot, with each new scene playing its
clip and every other shot showing the paper figure or notebook animation it starts from.

### 5. Put the whole trailer together, with the music

```bash
nice -n 19 python tools/trailer/cut.py --jobs 4        # frames, music, then build/magnus_trailer.mp4
```

`cut.py` draws every frame of every shot fresh at 1920×1080 (2906 frames, 96.9 s at 30 fps)
into `build/cut/`, using as many processes as `--jobs`; an interrupted run resumes, since frames
already drawn are skipped (delete `build/cut/` after changing a scene).  Then it writes
`build/music.wav` (`music.py`, also runnable on its own) and muxes both into
`build/magnus_trailer.mp4`: H.264 (CRF 18) with AAC audio at 192 kb/s.

What each kind of shot shows:

- a shot with a `scene`: that scene at the shot's progress;
- a pillar ("Accurate.", "Fast."): its scene with the big word and line over it;
- the code moment: the code types in, then the opening's matter curve draws from it;
- cards (the burden line, the reveal, "Beyond the textbook.", the ending): drawn in `cut.py`;
- a shot with `still_from`: **a placeholder until that shot is remade**, the paper figure (or the
  notebook animation, playing) framed on the dark ground with a slow push in, tagged
  "placeholder: to be remade".

Words go on a caption band at the bottom, several phrases one after another; cuts are hard, on
the beat; the picture fades in from black and out to black at the ends.

This step needs an ffmpeg with `libx264`, `aac` and PNG input, which the ffmpeg bundled with
Playwright lacks.  One that has them comes with `pip install imageio-ffmpeg`; point `$FFMPEG`
at `python -c "import imageio_ffmpeg; print(imageio_ffmpeg.get_ffmpeg_exe())"`.

### 6. Clean up

Everything generated lives in `tools/trailer/build/`, which git ignores; delete it to start
over.  At 1280×720 the frames average about 45 KB, so all eight scenes (1143 frames) take
about 50 MB; at 1920×1080 a frame is about 1.7 times larger.

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
