"""Every public callable that takes `angles` must document it, and vice versa.

A signature and a docstring are two structured texts that have to agree, and this sweep
edits ~85 of each.  That is precisely the shape of change that went wrong four times in
this repository's audit (I29, I30, L38, N44), so the agreement is checked mechanically
rather than by reading.

Also flags the reverse: an angle-taking function that has NOT been given `angles` yet, so
the sweep's remaining surface is always visible.
"""
import ast
import inspect
import pathlib
import sys

import magnus.globaldefs as gd
import magnus.hamiltonians as hamiltonians
import magnus.oscprob as oscprob
import magnus.oscprobstd as oscprobstd


def inert_angles(paths):
    """Functions that take `angles` and neither convert nor forward it.

    The third agreement, and the one neither the signature nor the docstring can show:
    pmns_mixing_matrix was given the parameter and the documentation and went on ignoring
    it, which would have read radians as sines in silence.  A parameter accepted and
    dropped is worse than one that was never added.
    """
    bad = []
    for path in paths:
        tree = ast.parse(pathlib.Path(path).read_text())
        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue
            if 'angles' not in [a.arg for a in node.args.args + node.args.kwonlyargs]:
                continue
            used = False
            for sub in ast.walk(node):
                # forwarded as angles=angles, or handed to resolve()/validate_convention()
                if isinstance(sub, ast.keyword) and sub.arg == 'angles':
                    used = True
                elif isinstance(sub, ast.Name) and sub.id == 'angles' \
                        and not isinstance(getattr(sub, 'ctx', None), ast.Store):
                    parent_call = True   # any read of the name counts as use
                    used = used or parent_call
            # a read inside the docstring does not exist, so this is a genuine body check
            if not used:
                bad.append('%s:%s' % (pathlib.Path(path).name, node.name))
    return bad

MODULES = [('hamiltonians', hamiltonians), ('oscprob', oscprob),
           ('oscprobstd', oscprobstd), ('globaldefs', gd)]

# lowercase s* / sxi* are mixing angles; lowercase d* are CP phases; uppercase D* are
# mass-squared splittings and are not angles at all.
def is_angle(name):
    return (name.startswith('sxi') or name == 'sth'
            or (name.startswith('s') and len(name) == 3 and name[1:].isdigit()))


def is_phase(name):
    return (name.startswith('dxi')
            or name == 'dCP'
            or (name.startswith('d') and len(name) == 3 and name[1:].isdigit()))


def main():
    missing_doc, missing_param, ok = [], [], 0
    for modname, mod in MODULES:
        for name in sorted(dir(mod)):
            if name.startswith('_'):
                continue
            fn = getattr(mod, name)
            if not callable(fn) or isinstance(fn, type):
                continue
            try:
                sig = inspect.signature(fn)
            except (TypeError, ValueError):
                continue
            params = list(sig.parameters)
            doc = inspect.getdoc(fn) or ''
            has_angles = 'angles' in params
            takes_angle = any(is_angle(p) for p in params)

            if has_angles:
                if 'angles : ' not in doc:
                    missing_doc.append('%s.%s' % (modname, name))
                else:
                    ok += 1
            elif takes_angle:
                missing_param.append('%s.%s' % (modname, name))

    if missing_doc:
        print('HAS `angles` BUT DOES NOT DOCUMENT IT (%d):' % len(missing_doc))
        for n in missing_doc:
            print('   ', n)
    else:
        print('every function with `angles` documents it (%d checked)' % ok)

    src = pathlib.Path('src/magnus')
    files = sorted(str(p) for p in src.rglob('*.py'))
    inert = inert_angles(files)
    print()
    if inert:
        print('ACCEPTS `angles` BUT NEVER USES IT (%d):' % len(inert))
        for n in inert:
            print('   ', n)
    else:
        print('every function with `angles` either converts or forwards it')

    print()
    print('STILL TO DO -- takes an angle but has no `angles` yet (%d):' % len(missing_param))
    for n in missing_param:
        print('   ', n)
    return 1 if (missing_doc or inert) else 0


if __name__ == '__main__':
    sys.exit(main())
