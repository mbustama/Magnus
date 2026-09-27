# The Magνs trailer

Everything needed to reproduce the trailer planned in issue #99, from the numbers to the
finished `magnus_trailer.mp4` with its music.  The trailer is a short advertisement for the
package, built from scenes computed with Magνs itself.  Every plot it shows is drawn here, in
the trailer's own look; no figure source in the repository is edited.

## Quick start

From the repository root:

```bash
pip install -e .                      # the package, NumPy, SciPy, Matplotlib (with Pillow)
pip install imageio-ffmpeg pymupdf    # an ffmpeg with H.264 and AAC; PDF thumbnails for the storyboard

nice -n 19 python tools/trailer/data.py                 # 1. every number the scenes show (a few minutes)
nice -n 19 python tools/trailer/cut.py --jobs 4         # 2. every frame, the music, the MP4
```

The result is `tools/trailer/build/magnus_trailer.mp4`: 1920×1080, 30 fps, H.264 video with
AAC stereo audio, 1:55 long (184 beats at 96 bpm).  The first run needs internet access
once, to download the two fonts (see "Fonts").

Everything the scripts write goes to `tools/trailer/build/`, which git ignores.

## The files

| File | What it is |
|---|---|
| `trailer.json` | The script and timeline: one entry per shot, in beats at 96 bpm, with its words (shown in the title band), what moves, the scene that draws it, and where its numbers come from.  Also the decisions still open. |
| `common.py` | Paths, the look (colors, fonts, easing), the physics setups of the new scenes, and an independent reference integrator (no Magnus code) that checks the adiabatic scene. |
| `data.py` | Computes every number the scenes show, with Magnus, and prints the checks. |
| `scenes.py` | The scenes made for the trailer: the two journeys, the burden, the Hamiltonian, how `strategy='auto'` picks a solver, the adiabatic switch, the flavors, the slab ladder, many energies in one call. |
| `paper_scenes.py` | The paper's figures remade as trailer scenes: the Earth oscillogram and its sterile-neutrino versions, Fermilab to four sites, CP violation, new physics, the flavor triangle, the Sun imaged, a buried body, geoneutrinos, a stellar jet, a long-range force. |
| `render.py` | Renders any scene alone, as a still, as PNG frames and as a clip.  Also finds ffmpeg for the other scripts. |
| `cut.py` | Assembles the whole trailer, shot by shot, as 1920×1080 frames, writes the music and muxes both into the MP4. |
| `music.py` | The music: an original piece synthesized from sine waves and noise. |
| `storyboard/` | `build.py` and `template.html`: a review page with one panel per shot. |

## Step by step

### 1. What you need

- **The package and its dependencies**: `pip install -e .` from the repository root.  The scenes
  were made with Matplotlib 3.10.8.
- **An ffmpeg with `libx264`, `aac` and PNG input**, for the MP4.  `pip install imageio-ffmpeg`
  provides one, and the scripts find it on their own.  To use another, set `$FFMPEG` to its
  path.  Without either, the scripts fall back to an `ffmpeg` on the PATH, and then to the one
  bundled with Playwright, which can write only silent WebM review clips: `cut.py` needs audio,
  so it needs one of the first two.
- **PyMuPDF** (`pip install pymupdf`), only for the storyboard's thumbnails of paper figures.
- **The paper's figure cache**, `notebooks/paper_figure_cache.json`, which is in the repository.
  `data.py` reads the paper's own arrays from it and never writes to it.

`nice -n 19` in the commands keeps a long render out of the way of other work.

### 2. Compute the data

```bash
nice -n 19 python tools/trailer/data.py                   # every step
nice -n 19 python tools/trailer/data.py paper adiabatic   # or only some
```

| Step | Writes | What it computes |
|---|---|---|
| `opening` | `build/opening.npz` | P(νμ→νe) at 3 GeV over 0–10,000 km, in vacuum and through an illustrative varying density |
| `code` | nothing | Checks that the code moment's user Hamiltonian runs without a warning and matches the shipped wrapper |
| `adiabatic` | `build/adiabatic.npz` | The adiabatic scene: 240 calls under `strategy='auto'` (about a minute), the levels, the resonances and the window Magnus patches, checked against the independent reference |
| `diagrams` | `build/diagrams.json` | The slab ladder at one point, a 2000-energy Earth spectrum, the mixing matrices for 2 to 5 flavors |
| `paper` | `build/paper.npz`, `build/land.json` | The remade paper scenes (seconds), and the continents of the paper's globes (notebook 28's `LAND` polygons, parsed from the notebook file): the oscillograms, the Sun, the jet, the geoneutrino curves and the long-range sweep from the paper's cache; Fermilab to four sites, CP violation, new physics, the flavor triangle and the buried body computed here with notebook 28's settings |

Each step prints a check line.  They should read:

```
opening: max P vacuum 0.108, matter 0.233
code moment: max |P - wrapper| = 3.7e-08, warnings: none
adiabatic: ... engines {'magnus': 15, 'hybrid': 225}, warnings none, max |P - reference| at 10 points = 5.2e-05, windows [(44993.6, 45005.1)]
diagrams: ladder [(2, 0.01118), ..., (21, 0.10203)], converged 0.10202; ...; |U|^2 rows sum to 1: True
paper: oscillograms (170, 200); Fermilab chords [763, 1284, 6724, 11632] km; cavity |dP| max 0.154; flavor triangle 122 + 61 points
```

### 3. Put the whole trailer together, with the music

```bash
nice -n 19 python tools/trailer/cut.py --jobs 4                 # frames, music, MP4
nice -n 19 python tools/trailer/cut.py --jobs 4 --frames-only   # frames only
nice -n 19 python tools/trailer/cut.py --jobs 4 --width 1280    # a smaller, faster review cut
```

`cut.py` works in three stages:

1. **Frames.** It draws every frame of every shot in `trailer.json` fresh, at 1920×1080 (3450
   frames at 30 fps), into `build/cut/f00000.png ...`, on as many processes as `--jobs`.  A run
   that is interrupted resumes where it stopped, since frames already drawn are skipped:
   **delete `build/cut/` after changing a scene or the script**, or old frames stay.
2. **Music.** It writes `build/music.wav` exactly as long as the picture (see "The music").
3. **Encoding.** It muxes the frames and the music into `build/magnus_trailer.mp4`:
   ```
   ffmpeg -framerate 30 -i build/cut/f%05d.png -i build/music.wav \
          -c:v libx264 -crf 18 -preset slow -pix_fmt yuv420p \
          -c:a aac -b:a 192k -movflags +faststart -shortest build/magnus_trailer.mp4
   ```

How each kind of shot is drawn:

- **A scene** (`"scene"` in the shot) is drawn at the shot's progress, in a box below the title
  band, with the shot's words in the band at the top of the frame (clear of a video player's
  controls), one phrase after another.  The last 12 % of a shot holds the finished picture.
  Each scene's finished picture is measured once (the rows of pixels that differ from the
  background) and raised so that its top sits just under the band; every frame of the shot uses
  that same offset, so nothing drifts from frame to frame.
- **A pillar** ("Accurate.", "Fast.") puts its big word and its line in the title band.
- **The code moment** types its code, then draws the opening's matter curve from it.
- **"Flexible."** and the **cards** (the question, the reveal, "From textbook to frontier.", the
  ending) are drawn in `cut.py` itself.

Cuts are hard, on the beat; the picture fades in from black and out to black at the ends.

### 4. Render one scene on its own (optional)

For working on a scene without assembling the whole cut:

```bash
nice -n 19 python tools/trailer/render.py switch earth_osc          # these scenes, 1280x720
nice -n 19 python tools/trailer/render.py --stills-only             # the finished frame of every scene
nice -n 19 python tools/trailer/render.py --width 1920 sun          # at full size
```

For each scene: `build/stills/<scene>.png` (the finished frame), `build/frames/<scene>/f*.png`
(every frame) and `build/clips/<scene>.mp4` (or `.webm` with Playwright's ffmpeg).  The
scenes are the keys of `SCENES` in `scenes.py` and `paper_scenes.py`; a scene's length here is
its default, and in the cut it takes its shot's length from `trailer.json`.

To change a scene, edit its function: it receives the figure, an axes spanning its box in a
16 × 9 coordinate system, and its progress `u` from 0 to 1.  `seg(u, a, b)` from `common.py`
gives the eased progress of a step that runs from `u = a` to `u = b`.  Axes it adds with
`fig.add_axes` take fractions of that box.

### 5. Review in the storyboard (optional)

```bash
python tools/trailer/storyboard/build.py        # -> build/storyboard.html
```

One self-contained page: the timeline, then a panel per shot with its words, what moves and
where its numbers come from, playing each scene's clip from `build/clips/` (run `render.py`
first) or showing its still.

### 6. Change the script

Edit `trailer.json`: a shot's `words` (the phrases of its title band, each starting with a
capital), its `beats` (its length at 96 bpm), or its `scene`.  Then delete `build/cut/` and run
`cut.py` again; the music follows the new timeline on its own, since `music.py` reads the acts
and their lengths from `trailer.json`.

### 7. Clean up

Delete `tools/trailer/build/` to start over.  The frames of the full cut are the bulk of it:
3450 frames at 1920×1080.

## The music

`music.py` synthesizes an original piece from sine waves, noise and a simple reverb: nothing is
sampled, so there is nothing to license.  It is deterministic (its noise uses a fixed seed), so
every run writes the same audio.  96 bpm in D minor, following each act's cue:

1. **Two journeys**: a quiet low drone, and a soft tone for each track.
2. **The problem**: a pulse enters.
3. **The reveal**: a low boom and the first full chord, then the groove.
4. **How it works**: a rising arpeggio over the groove.
5. **What you can do**: the full groove, building.
6. **From textbook to frontier**: the climax.  Its melody plays a computed probability: the
   Earth spectrum of the "Fast" diagram, P(νμ→νe) from 1 to 30 GeV, mapped onto the D-minor
   pentatonic scale.
7. **Accurate, fast, flexible**: three chord hits, resolving.
8. **Get it**: a final D-major chord with a long fade.

It is mastered to −1 dBFS.  Run it alone with `python tools/trailer/music.py`, which writes
`build/music.wav` for the length of the current script.

## Fonts

The trailer uses IBM Plex Mono and Unbounded, both under the SIL Open Font License.  On the
first run, `common.setup_matplotlib()` downloads them from the Google Fonts repository into
`build/fonts/`; after that nothing is fetched.  DejaVu, which ships with Matplotlib, fills in
the Greek letters Plex Mono lacks.

## What the scenes rest on

- **The two journeys**: P(νμ→νe) at 3 GeV over 0–10,000 km, in vacuum (`osc_prob_3nu_vacuum`) and
  through an illustrative density profile (`osc_prob_matter_std_potential`); NuFIT 6.1.
- **The code moment**: the user Hamiltonian shown on screen runs through
  `osc_prob_energy_baseline` without a warning and matches `osc_prob_matter_std_potential` to
  about 4e-8.
- **The adiabatic switch**: two flavors at 10 MeV crossing the MSW resonance slowly, through a
  shock, and slowly again.  One call per point under `strategy='auto'`: short paths are
  answered by the Magnus ladder and long ones by the hybrid, which finds the shock window
  itself; every checked point agrees with the independent reference to about 5e-5.  A sharper
  shock on a longer path is a case the hybrid gets wrong without a warning (issue #100); the
  scene stays out of that regime.
- **The diagrams**: the slab ladder is a real run at one point (fixed slab counts 2 to 21);
  the spectrum is a real 2000-energy `osc_prob_3nu_earth` call (no timing is claimed); the
  flavor bars are the package's own mixing matrices.
- **The globes** show the continents as the paper's do (notebook 28's `LAND` outlines, projected
  with the detector at the top of the limb: Kamioka, Fermilab, Gran Sasso).
- **The paper scenes**: the paper's own numbers (notebook 28), read from its cache where it
  cached them: the oscillograms (cell 42), the long-range force in the Sun (cell 62), the Sun
  imaged at five energies from 30 GeV to 3 TeV by a diffuse flux of neutrinos crossing it (cell 66, on
  its uniform impact-parameter grid),
  the jet (cell 69) and the geoneutrino curves (cells 86–88).  The rest is computed by
  `data.py paper` with the paper's settings: Fermilab to four sites (cell 44), new physics
  (cell 27), the flavor triangle (cell 84) and the buried body (cell 48, on 160 beam angles
  rather than the paper's 220).  The CP-violation ellipses are new: Fermilab to Homestake at
  2.5 GeV, both mass orderings, NuFIT 6.1.
