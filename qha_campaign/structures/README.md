# Starter structures for the QHA campaign

One primitive-cell `POSCAR_<phase>` per phase. The 4 here are **starter geometries** built to the
ambient V0 in `../../core/eos_params.json` and the standard prototype for each phase. They are
*starting points*: the static pre-flight relaxes ions + cell shape at fixed volume (ISIF=4), so
exact c/a self-corrects. The hcp-Fe validation gate (vs Dewaele 2006) catches any gross error.

| file | phase | prototype | atoms (prim) | matches manifest supercell |
|---|---|---|---|---|
| `POSCAR_hcp_Fe`  | hcp ε-Fe   | P6_3/mmc (A3)   | 2 | 4×4×3 → 96 |
| `POSCAR_B2_FeSi` | B2 FeSi    | Pm-3m (CsCl)    | 2 | 4×4×4 → 128 |
| `POSCAR_B1_FeO`  | B1 FeO     | Fm-3m (rocksalt)| 2 | 3×3×3 → 54 |
| `POSCAR_B8_FeO`  | B8 FeO     | P6_3/mmc (NiAs) | 4 | 3×3×2 → 72 |

## TWO phases you must supply (too complex to hand-write reliably)
- **Fe3S** — tetragonal `I-4` (Fe3P-type), 16 atoms. Source: Seagle et al. 2006 / Kamada 2014.
- **Fe7C3** — orthorhombic (Eckstrom-Adcock / `Pnma`) or hexagonal `P6_3mc`, 40 atoms.
  Source: Mookherjee/Nakajima 2011; the manuscript's highest-uncertainty phase.

Fetch a relaxed cell for each (Materials Project by formula, or a published CIF → POSCAR via
`pymatgen`/ASE), drop it here as `POSCAR_Fe3S` / `POSCAR_Fe7C3`, and confirm the primitive atom
count matches `make_manifest.py PHASES[...]['atoms_prim']` (16 and 40). If you change the cell,
update that dict and regenerate the manifest.

## Provenance note
These POSCARs are *own-derived starting structures* (prototype + sourced V0), not ingested from a
database, consistent with the repo's literature-only / own-derived discipline. The QHA F(V,T) that
comes out of them replaces the imported EoS block as `own_derived_parameter`.
