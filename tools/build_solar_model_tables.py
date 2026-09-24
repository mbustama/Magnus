# -*- coding: utf-8 -*-
r"""Build the solar-model tables that ship in ``src/magnus/data/solar_models/``.

Magnus ships twelve standard solar models, each trimmed to the three columns it reads: radius,
mass density (or its logarithm, where that is what the authors tabulate) and hydrogen mass
fraction.  The originals carry up to 47 columns and total about 11 MB, five times the rest of the
data the package needs; trimmed, they are about 0.7 MB.

Nothing is recomputed or reformatted.  Each value is the authors' own text, token for token, and
each table keeps the original file's header -- with its citation -- as comment lines, beneath a
provenance block naming the source, the date it was retrieved, the licence, and the SHA-256 of
the original file.  This script refuses to write anything from an original whose hash has
changed, so a table can always be traced to the exact file it came from.

The originals are not kept in the repository.  To rebuild, download them into one directory --
the URLs are in ``MODELS`` below -- and run, from the repository root::

    python tools/build_solar_model_tables.py DIRECTORY

Pass ``--check`` to verify instead that the committed tables are exactly what the originals give.
"""

import argparse
import hashlib
import pathlib
import sys
import zipfile

REPO = pathlib.Path(__file__).resolve().parents[1]
OUT = REPO / 'src' / 'magnus' / 'data' / 'solar_models'

BAHCALL_TERMS = ('"The files accessible here are available for general use. I would appreciate a note '
                 'at your convenience telling me how you are using them." (J. N. Bahcall, '
                 'https://www.sns.ias.edu/~jnb/SNdata/sndata.html)')
B16_TERMS = ('No terms stated on the authors\' page; cite Vinyoles et al. (2017).')
B23_TERMS = 'Creative Commons Attribution 4.0 International (CC-BY-4.0).'

B16_REF = ('N. Vinyoles, A. M. Serenelli, F. L. Villante, S. Basu, J. Bergström, M. C. Gonzalez-Garcia, '
           'M. Maltoni, C. Peña-Garay and N. Song, "A new generation of standard solar models", '
           'ApJ 835, 202 (2017), doi:10.3847/1538-4357/835/2/202, arXiv:1611.09867')
B16_PAGE = ('http://www.ice.csic.es/personal/aldos/Solar_Data.html (offline since at least 2026; the '
            'Internet Archive holds it)')
B23_REF = ('Y. Herrera and A. Serenelli, "Standard Solar Models B23 / SF-III", Zenodo, v1.2 (2023), '
           'doi:10.5281/zenodo.10822316')
B23_URL = 'https://zenodo.org/records/10822316/files/SolarModels_ICE2023_v1.2.zip'
B23_ZIP_SHA256 = '382037d76ed35c9b50a4094d3ada3c90284af025c3235b49cda9d932699f87d0'

# Composition references, exactly as the B23 README gives them.
B23_COMPOSITIONS = {
    'GS98': 'Grevesse & Sauval (1998), Space Sci. Rev., 85, 161',
    'AGSS09': 'Asplund et al. (2009), ARA&A, 47, 481',
    'C11': 'Caffau et al. (2011), Sol. Phys., 268, 255',
    'AAG21': 'Asplund et al. (2021), A&A 653, A141',
    'MB22m': 'Magg et al. (2022), A&A 661, A140 (meteoritic)',
    'MB22p': 'Magg et al. (2022), A&A 661, A140 (photospheric)',
}


def _bahcall(name, original, url, reference, retrieved, sha256):
    return dict(name=name, original=original, url=url, reference=reference, retrieved=retrieved,
                sha256=sha256, terms=BAHCALL_TERMS, layout='bahcall')


MODELS = [
    _bahcall('BP2000', 'bp2000stdmodel.dat',
             'https://www.sns.ias.edu/~jnb/SNdata/Export/BP2000/bp2000stdmodel.dat',
             'J. N. Bahcall, M. H. Pinsonneault and S. Basu, "Solar models: current epoch and time '
             'dependences, neutrinos, and helioseismological properties", ApJ 555, 990 (2001), '
             'doi:10.1086/321493, astro-ph/0010346',
             '2026-09-23', 'b84751e0f7008e6252df2e989b02827b2081d7652050296e0bb5c633a2ec9c14'),
    _bahcall('BP04', 'bp2004stdmodel.dat',
             'https://www.sns.ias.edu/~jnb/SNdata/Export/BP2004/bp2004stdmodel.dat',
             'J. N. Bahcall and M. H. Pinsonneault, "What do we (not) know theoretically about solar '
             'neutrino fluxes?", Phys. Rev. Lett. 92, 121301 (2004), '
             'doi:10.1103/PhysRevLett.92.121301, astro-ph/0402114',
             '2026-09-23', '91916dea000579d0982fba8ee2d4ef599e0e5603413de6975e5994a9d4ac245a'),
    _bahcall('BS05-OP', 'bs05op.dat',
             'https://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs05op.dat',
             'J. N. Bahcall, A. M. Serenelli and S. Basu, "New solar opacities, abundances, '
             'helioseismology, and neutrino fluxes", ApJ 621, L85 (2005), doi:10.1086/428929, '
             'astro-ph/0412440',
             '2026-09-23', '61c1600f4a221d8ed4b553349a3cd8b46fc0d893141d39baa8c63f9570237a2a'),
    _bahcall('BS05-AGS-OP', 'bs05_agsop.dat',
             'https://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs05_agsop.dat',
             'J. N. Bahcall, A. M. Serenelli and S. Basu, "New solar opacities, abundances, '
             'helioseismology, and neutrino fluxes", ApJ 621, L85 (2005), doi:10.1086/428929, '
             'astro-ph/0412440',
             '2026-08-04', '76da551f0a5c548f5f297269ed78937536c4238f2759ab101582658f9caf793e'),
    dict(name='B16-GS98', original='struct_b16_gs98.dat',
         url='https://web.archive.org/web/20220701130350/https://www.ice.csic.es/personal/aldos/'
             'Solar_Data_files/struct_b16_gs98.dat',
         reference=B16_REF, retrieved='2022-07-01 capture, from ' + B16_PAGE,
         sha256='ac9a1bdfd72699bd3d702d90b95ad853618b563475e072c9748a9e7e375f3f33', terms=B16_TERMS,
         layout='b16'),
    dict(name='B16-AGSS09met', original='struct_b16_agss09.dat',
         url='https://web.archive.org/web/20220701141642/https://www.ice.csic.es/personal/aldos/'
             'Solar_Data_files/struct_b16_agss09.dat',
         reference=B16_REF, retrieved='2022-07-01 capture, from ' + B16_PAGE,
         sha256='bd1add510bfe32f120ab6d858e25a0bbabfae11839bad0ce53d0da5b71ee7d2c', terms=B16_TERMS,
         layout='b16'),
] + [
    dict(name='B23-' + comp, original='SolarModels/struct+nu_SF3_%s.dat' % comp, url=B23_URL,
         reference=B23_REF + '; solar composition ' + ref, retrieved='2026-09-23',
         sha256=B23_ZIP_SHA256, terms=B23_TERMS, layout='b23', zipped='SolarModels_ICE2023_v1.2.zip')
    for comp, ref in B23_COMPOSITIONS.items()
]

# Which tokens of a data row are radius, density and hydrogen fraction, and what the density is.
LAYOUTS = {
    'bahcall': dict(columns=(1, 3, 6), density='rho_g_per_cm3', min_tokens=12),
    'b16': dict(columns=(1, 3, 6), density='rho_g_per_cm3', min_tokens=17),
    'b23': dict(columns=(0, 6, 18), density='log10_rho_g_per_cm3', min_tokens=47),
}


def _is_number(token):
    try:
        float(token)
    except ValueError:
        return False
    return True


def _read_original(model, directory):
    """The original file's text, after checking its hash (the zip's, for B23)."""
    directory = pathlib.Path(directory)
    if model.get('zipped'):
        path = directory / model['zipped']
        blob = path.read_bytes()
        if hashlib.sha256(blob).hexdigest() != model['sha256']:
            raise SystemExit('%s: SHA-256 does not match; refusing to build from it' % path)
        return zipfile.ZipFile(path).read(model['original']).decode('utf-8')
    path = directory / model['original']
    blob = path.read_bytes()
    if hashlib.sha256(blob).hexdigest() != model['sha256']:
        raise SystemExit('%s: SHA-256 does not match; refusing to build from it' % path)
    return blob.decode('utf-8')


def trimmed_table(model, directory):
    """The text of the trimmed table for ``model``."""
    layout = LAYOUTS[model['layout']]
    text = _read_original(model, directory).replace('\r\n', '\n').replace('\r', '\n')
    header, rows = [], []
    for line in text.split('\n'):
        tokens = line.split()
        if len(tokens) >= layout['min_tokens'] and all(_is_number(t) for t in tokens):
            rows.append([tokens[i] for i in layout['columns']])
        elif not rows:
            header.append(line.rstrip())            # everything before the table, verbatim
    if model['layout'] == 'b23':
        header = header[:1]                          # the column line; the package README is cited instead
    while header and not header[-1].strip():
        header.pop()
    out = [
        '# Magnus solar model %s' % model['name'],
        '#',
        '# Reference:  %s' % model['reference'],
        '# Source:     %s' % model['url'],
        '# Original:   %s' % model['original'],
        '# Retrieved:  %s' % model['retrieved'],
        '# SHA-256:    %s (of %s)' % (model['sha256'], model.get('zipped', model['original'])),
        '# Terms:      %s' % model['terms'],
        '# Modified:   columns extracted only -- radius [R_sun], %s, hydrogen mass fraction X --'
        % layout['density'],
        '#             each value copied from the original as written; nothing recomputed or rounded.',
        '#             Built by tools/build_solar_model_tables.py.',
        '# density_column: %s' % layout['density'],
        '#',
        '# The original header follows, as the authors wrote it:',
        '#',
    ]
    out += ['#   ' + h if h.strip() else '#' for h in header]
    out += ['#', '# r_over_r_sun  %s  x_hydrogen' % layout['density']]
    out += ['  '.join(r) for r in rows]
    return '\n'.join(out) + '\n'


def file_name(model):
    return model['name'].lower().replace('-', '_') + '.dat'


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('directory', help='directory holding the original files')
    ap.add_argument('--check', action='store_true',
                    help='verify the committed tables instead of writing them')
    args = ap.parse_args()
    stale = []
    OUT.mkdir(parents=True, exist_ok=True)
    for model in MODELS:
        text = trimmed_table(model, args.directory)
        path = OUT / file_name(model)
        if args.check:
            if not path.exists() or path.read_text(encoding='utf-8') != text:
                stale.append(path.name)
        else:
            path.write_text(text, encoding='utf-8')
            n = sum(1 for line in text.split('\n') if line and not line.startswith('#'))
            print('%-22s %5d rows  %8d bytes  -> %s' % (model['name'], n, len(text.encode()),
                                                        path.relative_to(REPO)))
    if args.check:
        print('stale: %s' % (', '.join(stale) or 'none'))
        return 1 if stale else 0
    return 0


if __name__ == '__main__':
    sys.exit(main())
