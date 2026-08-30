# The paper

The source of *Magνs: neutrino oscillation probabilities for any Hamiltonian,
any number of flavors, any density profile*, the Computer Physics
Communications article that documents this library.

It lives here so that the paper and the code it describes travel together: the
figures in it are produced by `notebooks/28_magnus_paper_figures.ipynb`, from
data computed as that notebook runs where the number is Magνs's own and from
the frozen `notebooks/external_*.json` datasets where it belongs to another
code. A claim in the text can therefore be traced to the measurement behind it
without leaving the repository.

Nothing here is installed. The wheel is built from `src/`, and the source
distribution from setuptools' default set; neither reaches `resources/`. It is
here for a clone, not for `pip`. `tests/test_file_tree.py` knows about these
files, so adding one to this folder means updating `TREE` there — see the note
at the end.

## Building it

```bash
cd resources/paper
pdflatex main.tex && bibtex main && pdflatex main.tex && pdflatex main.tex
```

Four passes: once to collect the citations, `bibtex` to resolve them, and twice
more for the cross-references to settle. **This differs from the sibling
project's README**, which still says there is no `bibtex` step — that was true
of an older revision of that paper, whose bibliography was an inlined
`thebibliography`. Its current `main.tex` carries no `\bibitem` at all and runs
`\bibliography{refs.bib}` live, and so does this one. A new citation goes into
`refs.bib` and is picked up on the next build; nothing is numbered by hand.

`refs.bib` is the NuOscProbExact bibliography with the Magνs-specific entries
appended below a marked separator. Carrying the whole of it means that a
reference cited in both papers has identical metadata in both, and BibTeX emits
only what is cited, so the unused majority costs nothing.

Everything the preamble loads is in a normal TeX Live install. `elsarticle.cls`
and `elsarticle-num.bst` are bundled anyway, so that the folder compiles on a
machine whose TeX Live was installed without the Elsevier bundle.

## Line numbers

On, for a referee, through a switch in the preamble:

```latex
\linenostrue         % \linenosfalse for the camera-ready
```

`elsarticle`'s own `review` option will *not* do this — it only sets preprint
mode and 1.5 spacing, and discards the `5p` two-column layout — so `lineno` is
loaded directly. `switch` puts the numbers on the outer edge of each column,
which two-column needs: without it the right column's numbers land in the
gutter. `mathlines` numbers display equations. Floats are never numbered, so
the figures, tables and listings are skipped.

## Revisions

**Write ordinary LaTeX in `main.tex`.** Nothing in it should record what changed
since a previous version; that is worked out mechanically with `latexdiff`
against a frozen baseline of the submitted version. There is no baseline yet —
this is a first submission — so freeze `main.tex` as `baseline_cpc_v1.tex` when
it is submitted, and diff against that afterwards.

The sibling project's `make_versions.py` records five things that had to be
worked around to get `latexdiff` through a document of this shape, each of
which failed in a way that did not point at its own cause. Read it before
setting the diff up here rather than rediscovering them:

- the bibliography must be held out of the diff entirely, since a marked-up URL
  inside a `\bibitem` sends `hyperref` into a recursion that exhausts TeX's
  input stack;
- `\texorpdfstring` has to come out of the two mark-up macros;
- `--graphics-markup=none`, since the default boxing of changed figures
  recurses inside a caption in a `figure*`;
- `--disable-citation-markup`, which is cosmetic rather than fatal.

## The figures

All ten are in `figs/`, found through a single `\graphicspath{{figs/}}` so
that no `\includegraphics` names the folder. All ten are produced by
`notebooks/28_magnus_paper_figures.ipynb`. To rebuild them:

```bash
python notebooks/make_notebooks.py --only 28
```

which writes the notebook, executes it, and stores its outputs. About ninety
seconds.
The notebook writes into `resources/paper/figs/` by default; set
`MAGNUS_PAPER_FIGDIR` to send them elsewhere.

Note that this rewrites every PDF whether or not anything changed, and the
bytes differ between runs, so `git status` will show them as modified
afterwards even when the figures are identical. Commit them only when a figure
actually changed.

The comparison figures draw their external-code numbers from data frozen in
`notebooks/`, so none of the other codes has to be installed to redraw them:

| file | what it holds |
|---|---|
| `external_profile_benchmarks.json` | the smooth-profile plane, 2 to 5 flavors |
| `external_prem_speed_accuracy.json` | the two Earth planes, five external codes |
| `external_shock_benchmarks.json` | the supernova shock at two front widths |
| `external_solar_nusquids.json` | nuSQuIDS on the same solar model file |
| `external_speed_accuracy.json` | the six-code constant-density plane |

One further dataset in `notebooks/` is not another code's but our own, and is
frozen for a different reason. The right panel of Fig. 2 measures deviations
down to $10^{-15}$, which no double-precision reference can certify, so it is
measured against a midpoint slab product carried at fifty decimal digits in
`mpmath` and Richardson-extrapolated three times. That reference is converged to
$3 \cdot 10^{-20}$ — nothing in the panel is the reference rather than the
method — and it costs about two and a half minutes to build, which is not worth
paying on every rebuild of the figure. It is therefore stored in
`mpmath_phase_reference.json` under a fingerprint of the configuration that
produced it: nine samples of the Hamiltonian along the trajectory, the baseline,
the working precision and the slab counts. Change any of those and the notebook
recomputes and rewrites the file by itself; change none of them and it reads it.
Nothing has to be invalidated by hand.

Magνs's own numbers on those planes are computed live rather than frozen,
because a timing is only comparable if it was taken on the machine the other
timings came from — and because the conventions that make the comparison mean
anything (the matched matter potential above all) are then visible in the
notebook rather than buried in a generator that ran once.

## Timings

Every timing quoted in the paper was measured on one machine: an Intel Core
i5-1334U, ten cores and twelve threads to 4.6 GHz, 16 GB of memory, Ubuntu
24.04.4 LTS on kernel 7.0.0-29-generic, with Python 3.12.7, NumPy 1.26.4,
SciPy 1.15.3 and Numba 0.60.0. Absolute numbers will differ elsewhere; the
ratios between codes and between settings, which is what the figures are about,
are far more stable than that. Re-measure on an idle machine, or interleave
with a control workload the change cannot touch — the package's own timing
harness does the latter, and it is what makes its ratios readable at all.

## Adding a file here

`tests/test_file_tree.py` requires the tracked files to match the documented
tree exactly, so a new file fails the suite until the table follows it:

```bash
git add resources/paper/<newfile>
python3 tests/test_file_tree.py --write
```

`figs/` is collapsed to a single entry in that table, as `img/gallery/` and
`docs/dev/` are, so a new *figure* needs no `TREE` update — but it does need a
`git add -f`. The repository gitignores `*.pdf` globally, and the ten figures are
tracked past that ignore on purpose, so that a clone compiles the paper without
first executing a notebook:

```bash
git add -f resources/paper/figs/<newfigure>.pdf
```
