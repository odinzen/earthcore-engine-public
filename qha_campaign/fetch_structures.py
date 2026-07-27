"""
Pull a structure from Materials Project and write it as a POSCAR for the campaign.
Needs `mp-api` in the env and MP_API_KEY set (free key from materialsproject.org -> API).

  # 1. list the polymorphs for a formula, pick the right prototype (space group + atom count):
  python fetch_structures.py --list Fe3S
  python fetch_structures.py --list Fe7C3

  # 2. fetch the chosen one to a POSCAR:
  python fetch_structures.py --get mp-XXXXX --out structures/POSCAR_Fe3S

Pick by PHYSICS, not just lowest energy: Fe3S wants the Fe3P-type tetragonal I-4 (Seagle 2006,
Kamada 2014 -- the core sulfide); Fe7C3 wants the orthorhombic (Pnma) or hexagonal P6_3mc cell
(Mookherjee 2011). The list shows space group + n_sites so you can match.
"""
import argparse
import os
import sys

KEY = os.environ.get("MP_API_KEY")


def client():
    if not KEY:
        sys.exit("MP_API_KEY not set (get a free key at materialsproject.org -> API)")
    try:
        from mp_api.client import MPRester
    except ImportError:
        sys.exit("mp-api not installed:  mamba install -c conda-forge mp-api -y")
    return MPRester(KEY)


def do_list(formula):
    with client() as mpr:
        docs = mpr.materials.summary.search(
            formula=formula,
            fields=["material_id", "formula_pretty", "symmetry", "nsites",
                    "energy_above_hull", "volume"])
    if not docs:
        sys.exit(f"no Materials Project entries for {formula}")
    print(f"\n{formula}: {len(docs)} polymorph(s) (pick by space group + n_sites)")
    print(f"  {'material_id':14} {'spacegroup':14} {'nsites':>6} {'eV/at_hull':>11} {'A3/atom':>9}")
    for d in sorted(docs, key=lambda x: (x.energy_above_hull if x.energy_above_hull is not None else 9)):
        sg = getattr(d.symmetry, "symbol", "?")
        eah = d.energy_above_hull if d.energy_above_hull is not None else float("nan")
        print(f"  {d.material_id!s:14} {sg:14} {d.nsites:6d} {eah:11.3f} {d.volume/d.nsites:9.2f}")


def do_get(mpid, out):
    from pymatgen.io.vasp import Poscar
    with client() as mpr:
        struct = mpr.get_structure_by_material_id(mpid)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    Poscar(struct).write_file(out)
    print(f"wrote {out}: {struct.composition.reduced_formula}, {len(struct)} atoms, "
          f"{struct.get_space_group_info()}")
    print("  -> confirm n_atoms matches make_manifest.py PHASES; regenerate the manifest if changed.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", metavar="FORMULA")
    ap.add_argument("--get", metavar="MPID")
    ap.add_argument("--out", default="structures/POSCAR_out")
    a = ap.parse_args()
    if a.list:
        do_list(a.list)
    elif a.get:
        do_get(a.get, a.out)
    else:
        ap.error("use --list FORMULA or --get MPID --out PATH")


if __name__ == "__main__":
    main()
