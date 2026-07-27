"""Dump the QHA inputs (static energy + F_vib(T) per volume) to a small JSON per phase,
so the EoS fit and its validation gates can run and be checked off-cluster. phonopy does the
heavy lifting here (on Sol); everything downstream needs only numpy.

    python dump_qha.py --phase hcp_Fe

Writes qha_<phase>.json with the temperature grid and, per complete volume, the cell volume,
the static total energy, and F_vib(T). Incomplete volumes are skipped with a note, matching
fit_qha. Reuses fit_qha for the phonon free energy so the numbers are identical to the inventory.
"""
import argparse
import json
import os

import numpy as np

import fit_qha as F

RY_TO_EV = 13.605693122994


def parse_energy(pwout):
    """Total energy (eV) from the last '!' line of a QE output (the converged value)."""
    e = None
    for line in open(pwout):
        if line.lstrip().startswith("!") and "total energy" in line:
            e = float(line.split("=")[1].split()[0])          # Ry
    if e is None:
        raise SystemExit(f"{pwout}: no '! total energy' line")
    return e * RY_TO_EV


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", required=True)
    a = ap.parse_args()

    temps = np.arange(F.T_MIN, F.T_MAX + F.T_STEP, F.T_STEP)   # 300..7000 K, matches fit_qha
    vols = []
    for v in F.volumes_for(a.phase):
        static = os.path.join(F.HERE, "runs", a.phase, f"v{v:02d}", "static", "pw.out")
        try:
            V, fvib = F.free_energy(a.phase, v)                # V (Ang^3), F_vib(T) eV/cell
            E = parse_energy(static)
        except SystemExit as e:
            print(f"[dump] {a.phase} v{v:02d}: SKIP -- {e}")
            continue
        if len(fvib) != len(temps):
            raise SystemExit(f"{a.phase} v{v:02d}: F_vib has {len(fvib)} points, "
                             f"expected {len(temps)} -- T-grid mismatch, tell Michael")
        vols.append({"vol_idx": v, "V_ang3": float(V), "E_static_eV": float(E),
                     "F_vib_eV": [float(x) for x in fvib]})
        print(f"[dump] {a.phase} v{v:02d}: V={V:.2f} A^3  E={E:.4f} eV  F_vib(300K)={fvib[0]:.4f} eV")

    out = os.path.join(F.HERE, f"qha_{a.phase}.json")
    with open(out, "w") as fh:
        json.dump({"phase": a.phase, "temperatures_K": [float(t) for t in temps], "volumes": vols},
                  fh, indent=1)
    print(f"[dump] wrote {out}: {len(vols)} complete volumes")


if __name__ == "__main__":
    main()
