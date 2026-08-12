"""Does every documented "Default: X" match the actual signature default?

A docstring that states a default is a claim, and this package has already shipped two
defaults that moved (the implicit NuFit release, the Earth's Y_e) without every mention
following.  This compares the two mechanically rather than by reading.
"""
import inspect
import re
import sys

import magnus.earth as earth
import magnus.globaldefs as gd
import magnus.hamiltonians as hams
import magnus.magnus as mg
import magnus.matter as matter
import magnus.oscprob as op
import magnus.oscprobstd as ops

MODULES = [hams, op, ops, gd, matter, earth, mg]

# "name : type, optional" then prose containing "Default: <value>."
ENTRY = re.compile(r'^\s*(\w+)\s*:\s*[^\n]*\n((?:\s{4,}.*\n?)*)', re.M)
# to end of line, then strip one trailing period -- splitting on the first '.' turned
# 0.5 into 0 and './out.log' into a lone quote, which is 150 false positives
DEFAULT = re.compile(r'Default:\s*([^\n]+)')


def normalise(text):
    t = text.strip()
    if t.endswith('.'):
        t = t[:-1]
    t = t.strip().strip('`').strip("'").strip('"')
    return {'True': 'True', 'False': 'False', 'None': 'None'}.get(t, t)


def same(doc_val, real):
    d = normalise(doc_val)
    if d in ('True', 'False', 'None'):
        return d == repr(real)
    try:
        return abs(float(d) - float(real)) <= 1e-12*max(1.0, abs(float(real)))
    except (TypeError, ValueError):
        return d == str(real) or d.strip("'\"") == str(real)


def main():
    bad = []
    seen = set()
    for mod in MODULES:
        for name in sorted(dir(mod)):
            if name.startswith('_'):
                continue
            fn = getattr(mod, name)
            if not callable(fn) or isinstance(fn, type):
                continue
            key = getattr(fn, '__module__', '') + '.' + name
            if key in seen:
                continue
            seen.add(key)
            try:
                sig = inspect.signature(fn)
            except (TypeError, ValueError):
                continue
            doc = inspect.getdoc(fn) or ''
            if 'Parameters' not in doc:
                continue
            for m in ENTRY.finditer(doc):
                pname, body = m.group(1), m.group(2)
                if pname not in sig.parameters:
                    continue
                real = sig.parameters[pname].default
                if real is inspect.Parameter.empty:
                    continue
                dm = DEFAULT.search(body)
                if not dm:
                    continue
                doc_val = dm.group(1).strip()
                if real is None and ('None' in doc_val or ':data:' in doc_val
                                     or 'layered' in doc_val or 'every ' in doc_val):
                    continue      # prose describing what None resolves to, not a literal
                if not same(doc_val, real):
                    bad.append((key, pname, dm.group(1).strip(), repr(real)))

    if bad:
        print('DOCUMENTED DEFAULT DISAGREES WITH THE SIGNATURE (%d):' % len(bad))
        for k, pn, doc_v, real in bad:
            print('  %-52s %-22s doc=%-18s real=%s' % (k, pn, doc_v, real))
    else:
        print('every documented default matches its signature')
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
