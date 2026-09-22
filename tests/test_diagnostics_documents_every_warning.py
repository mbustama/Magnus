# -*- coding: utf-8 -*-
r"""Every warning the package can raise appears in the diagnostics catalogue.

``docs/source/diagnostics.rst`` opens by promising "what every warning means", and
its catalogue had nine of the fourteen classes.  Four of the five missing ones --
``BaselineUnitWarning``, ``CrossCheckInconclusiveWarning``,
``PseudoDiracSplittingWarning`` and ``SterileMatterCompositionWarning`` -- appeared
nowhere in ``docs/source`` at all, and the last of those carries the caveat every
3+1 user needs about the sterile states' entry in the matter projector.

A new warning class is easy to add and easy to forget to document, which is how
that happened; this fails when it does.  Notebook 20's table is the other place
that lists them all, and is the natural source for a new row's wording.
"""

import importlib
import inspect
import pkgutil
import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DIAGNOSTICS = REPO / 'docs' / 'source' / 'diagnostics.rst'


def _warning_classes():
    """Every Warning subclass magnus defines, including the nested subpackage's."""
    import magnus
    found = {}
    modules = ['magnus'] + ['magnus.' + m.name
                            for m in pkgutil.walk_packages(magnus.__path__)]
    for name in list(modules):
        # walk_packages only reaches the top level; hamiltonians is a subpackage and
        # hides PseudoDiracSplittingWarning, which a flat scan misses.
        try:
            mod = importlib.import_module(name)
        except Exception:
            continue
        modules.extend(name + '.' + m.name
                       for m in pkgutil.walk_packages(getattr(mod, '__path__', [])))
        for obj in vars(mod).values():
            if (inspect.isclass(obj) and issubclass(obj, Warning)
                    and obj.__module__.startswith('magnus')):
                found[obj.__name__] = obj.__module__
    return found


def test_every_warning_class_is_in_the_catalogue():
    defined = _warning_classes()
    assert len(defined) >= 14, "found only %d warning classes: %s" % (
        len(defined), sorted(defined))

    text = DIAGNOSTICS.read_text(encoding='utf-8')
    documented = set(re.findall(r'class:`magnus\.[\w.]*?\.(\w+Warning)`', text))

    missing = sorted(set(defined) - documented)
    assert not missing, (
        "%s promises 'what every warning means' and does not mention: %s.\n"
        "Add a row to its catalogue table; notebook 20 lists all of them with "
        "the wording to reuse." % (DIAGNOSTICS.relative_to(REPO), ', '.join(missing)))

    stale = sorted(documented - set(defined))
    assert not stale, (
        "%s documents warnings the package no longer defines: %s"
        % (DIAGNOSTICS.relative_to(REPO), ', '.join(stale)))
