# -*- coding: utf-8 -*-
r"""Every environment variable the documentation tells CI to set, CI actually sets.

`MAGNUS_PAPER_CACHE_ONLY` was written together with `notebooks/paper_figure_cache.json`,
and `resources/paper/README.md` has instructed that continuous integration set it ever
since.  No workflow was ever told.  `git log -S` over every branch finds the name nowhere
under `.github/`, so the rule held only on machines where somebody exported it by hand.

The cost of that gap is not mainly the hour a recompute takes.  It is that a cache miss
**passed**: the notebook recomputed inputs that had drifted and redrew the figure from
them, and the gate reported green.  A stale oracle that does not say so is the failure
this very branch exists to fix, and this was the same shape of it one level up.

So the instruction is now executable.  A sentence in the documentation that names a
`MAGNUS_*` variable and the phrase "continuous integration" is read as a requirement, and
the workflows have to satisfy it.  Writing the rule down again is not what stops this
recurring; checking it is.
"""

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parent.parent
WORKFLOWS = ROOT/'.github'/'workflows'
# Where an instruction to CI may be written.  Kept short deliberately: a prose file that
# merely mentions a variable is not giving an order, and a test that guesses otherwise
# fails for reasons nobody can act on.
INSTRUCTION_FILES = ('resources/paper/README.md',)


def uncommented(path):
    """The workflow's lines with whole-line YAML comments removed.

    A variable set on a commented-out line is not set.  That is the one way this test
    could pass while the thing it checks is untrue, so it is the one thing it strips.
    """
    return [line for line in path.read_text().splitlines()
            if not line.lstrip().startswith('#')]


def required_variables():
    """Every MAGNUS_* variable named in a sentence that also says continuous integration."""
    wanted = {}
    for name in INSTRUCTION_FILES:
        text = (ROOT/name).read_text()
        for sentence in re.split(r'(?<=[.!?])\s+', text):
            if 'continuous integration' not in sentence.lower():
                continue
            for var in re.findall(r'MAGNUS_[A-Z0-9_]+', sentence):
                wanted.setdefault(var, name)
    return wanted


def test_documentation_names_at_least_one_requirement():
    """The extraction above is not silently matching nothing.

    Without this, rewording the README into a form the regex misses would empty the rule
    and every other test here would pass by vacuity -- which is how the gap being closed
    went unnoticed in the first place.
    """
    assert required_variables(), (
        'no MAGNUS_* variable was found in a sentence about continuous integration in %s; '
        'either the instruction was removed or this test no longer recognises it'
        % ', '.join(INSTRUCTION_FILES))


def test_every_documented_ci_variable_is_set_in_a_workflow():
    """What the documentation orders, a workflow carries out."""
    lines = [line for path in sorted(WORKFLOWS.glob('*.yml'))
             for line in uncommented(path)]
    for var, source in sorted(required_variables().items()):
        assert any(re.search(r'\b%s\s*:\s*\S' % re.escape(var), line) for line in lines), (
            '%s instructs continuous integration to set %s, and no workflow in %s does. '
            'Either set it on the job that needs it, or stop saying so.'
            % (source, var, WORKFLOWS.relative_to(ROOT)))


def test_the_paper_figures_are_rendered_and_never_recomputed():
    """The specific rule, asserted where a reader of the workflow will look for it.

    Named separately from the general test because this one carries the reason: the
    figures' inputs are machine-specific timings and settled references, so a shared
    runner re-deriving them learns nothing and reports drift as agreement.
    """
    notebooks = WORKFLOWS/'notebooks.yml'
    assert any(re.search(r'\bMAGNUS_PAPER_CACHE_ONLY\s*:\s*\S', line)
               for line in uncommented(notebooks)), (
        'notebooks.yml does not set MAGNUS_PAPER_CACHE_ONLY, so a paper-figure cache miss '
        'on the runner will recompute quietly and report a drifted figure as a pass')
