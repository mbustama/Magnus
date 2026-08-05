# Third-party data in this directory

One file here is **not** ours. This note records where it came from, what its terms are, and
how to replace it, so that anyone redistributing this repository knows what they are
redistributing.

## `bs05_agsop.dat` — the BS2005-AGS,OP standard solar model

| | |
|---|---|
| **What** | Tabulated standard solar model: mass fraction, radius, temperature, density, pressure, luminosity fraction, and the mass fractions of ¹H, ⁴He, ³He, ¹²C, ¹⁴N and ¹⁶O, on 1284 radial zones from 0.0016 to 0.983 R☉ |
| **Source** | <https://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs05_agsop.dat> |
| **Retrieved** | 2026-08-04, 158 748 bytes, unmodified |
| **Citation** | J. N. Bahcall, A. M. Serenelli & S. Basu, *"New solar opacities, abundances, helioseismology, and neutrino fluxes"*, ApJ **621**, L85 (2005) — astro-ph/0412440 |
| **Stated terms** | The source page states no licence and no copyright restriction. Its only usage statement is a courtesy request, quoted verbatim below. |

> "The files accessible here are available for general use. I would appreciate a note at your
> convenience telling me how you are using them."

So redistribution is explicitly permitted. There is no licence conflict with this repository's
GPL-3.0: that licence governs the code here, and this file is third-party *data* carried
alongside it, not a derived work of it.

**The courtesy request is worth honouring** if you build on this. Note that John Bahcall died in
2005; the standard solar models are now maintained by Aldo Serenelli, and current tables live at
his pages rather than the IAS mirror above.

## Why the file is committed rather than downloaded

Reproducibility. Every diagnostic in this directory exists so that a number in
`../FINDINGS_ROBUSTNESS_PROGRAMME.md` can be re-derived years later, and a script that silently
depends on a URL staying alive is not reproducible — the IAS mirror already returns a 302 for
the `http://` form, and pages of this age disappear. 155 KB is a small price for a result that
still reproduces when the host does not.

The same reasoning is why the physical-profile population carries a `provenance` string per
family (`physical_profiles.py`): a finding is only as good as the profile it was measured on,
and "where did this shape come from" has to survive the session that produced it.

## Re-fetching it

```bash
curl -L -o bs05_agsop.dat https://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/bs05_agsop.dat
```

`physical_profiles.load_bs05()` parses it by requiring twelve float fields per line, which skips
the header, the column-heading line and the two trailing `Lsun=`/`Rsun=` lines. It computes the
electron number density as `n_e = rho * N_A * (1 + X) / 2` — the hydrogen mass fraction is read
from the table rather than assumed, because `Y_e` runs from 0.68 at the centre to 0.86 at the
surface and a fixed 0.5 would be wrong by up to 70 %.

## Everything else here

Every other file in this directory is original to this repository. The **shapes** of two
profile families are taken from the literature and are cited in `physical_profiles.py`'s
`References` block — Fogli, Lisi, Mirizzi & Montanino (PRD **68**, 033005) for the supernova
progenitor profile and shock parametrization, and Kneller & Kabadi (PRD **92**, 013009) for the
turbulence construction and the shock radii — but the implementations are ours, and no code or
data from those papers is included.
