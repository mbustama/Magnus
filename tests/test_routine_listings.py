# -*- coding: utf-8 -*-
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (C) 2026 Mauricio Bustamante
r"""Every module's ``Routine listings`` has to keep up with the module.

Each module docstring carries a ``Routine listings`` section naming what the module
provides, and it is the first thing a reader meets on the rendered API page.  Nothing
regenerated it, so it drifted silently: the pre-publish audit added
``earth.electron_fraction_func_prem``, ``earth.neutron_to_proton_ratio_from_electron_fraction``,
``matter.matter_potential_projector`` and ``adiabatic.oscillation_sampling`` and listed none
of them, and that was only noticed by reading the modules against their own docstrings months
later.

This is deliberately lenient about *format* -- the sections use two different bullet styles
already -- and strict about *coverage*: a public function defined in a module has to be
mentioned somewhere in that module's listing section.  Naming it is the point; how it is
punctuated is not.

Classes are not required.  The sections are titled "Routine listings" and the modules that
define warning classes describe them in prose instead, which is a defensible reading of
"routine".
"""

import importlib
import inspect
import pathlib
import re

import pytest


ROOT = pathlib.Path(__file__).resolve().parent.parent
SRC = ROOT/'src'

# The two modules that deliberately carry no listing, pinned so that a section cannot
# quietly disappear from a module that has one: dropping one here is a visible edit.
#
#   oscprob  -- ninety-odd public wrappers, generated in families; a hand-kept list of them
#               would be longer than the module docstring and stale within a release.
#   plotting -- its docstring is organised around the figure types, not the call names.
NO_LISTING = {'magnus.oscprob', 'magnus.plotting'}


def _modules():
    out = []
    for path in sorted(SRC.rglob('*.py')):
        if path.name == '__main__.py':
            continue
        out.append('.'.join(path.relative_to(SRC).with_suffix('').parts))
    return out


def _listing_section(doc):
    """The text from the ``Routine listings`` heading to the end of the docstring."""
    m = re.search(r'Routine listings\s*\n-+\n', doc)
    return doc[m.end():] if m else None


def _public_functions(mod):
    """Public functions *defined here*, so re-exports are the defining module's business."""
    return sorted(
        name for name, obj in vars(mod).items()
        if not name.startswith('_')
        and inspect.isfunction(obj)
        and getattr(obj, '__module__', None) == mod.__name__)


@pytest.mark.parametrize('modname', _modules())
def test_routine_listings_name_every_public_function(modname):
    mod = importlib.import_module(modname)
    doc = inspect.getdoc(mod) or ''
    section = _listing_section(doc)

    if modname in NO_LISTING:
        assert section is None, (
            '%s is listed in NO_LISTING but now has a Routine listings section; remove it '
            'from that set so the section is kept complete.' % modname)
        return

    assert section is not None, (
        '%s has no Routine listings section.  Add one, or add the module to NO_LISTING '
        'with a reason.' % modname)

    missing = [f for f in _public_functions(mod)
               if not re.search(r'\b%s\b' % re.escape(f), section)]
    assert not missing, (
        '%s defines public functions its Routine listings does not name: %s'
        % (modname, ', '.join(missing)))
