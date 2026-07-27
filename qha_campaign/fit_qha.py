"""
Post-processing: forces -> phonon free energy -> quasiharmonic F(V,T) -> own EoS.
Reads forces and volume DIRECTLY from the QE output (not via ASE's full parser, which asserts on
eigenvalues that verbosity=low omits). Uses phonopy's Python API. eV/Ang throughout.

    python fit_qha.py --phase B1_FeO

Per phase, for each volume: load disp/phonopy_disp.yaml, read supercell-*.out forces, build force
constants, thermal properties -> F_vib(V,T). (EoS fit -> eos_params_qha.json still TODO.)
"""
import argparse
import csv
import os

import numpy as np
import phonopy

HERE = os.path.dirname(os.path.abspath(__file__))
T_MIN, T_MAX, T_STEP = 300, 7000, 100
RY_BOHR_TO_EV_ANG = 25.71104309
BOHR3_TO_ANG3 = 0.529177210903 ** 3
KJMOL_TO_EV = 1000.0 / 96485.33212


def volumes_for(phase):
    vs = set()
    with open(os.path.join(HERE, "job_manifest.csv")) as f:
        for r in csv.DictReader(f):
            if r["phase"] == phase:
                vs.add(int(r["vol_idx"]))
    return sorted(vs)


def parse_forces(pwout):
    """Forces (eV/Ang) from the last 'Forces acting on atoms' block of a QE output."""
    if not os.path.isfile(pwout):
        raise SystemExit(f"{pwout}: missing -- phonon job not finished")
    lines = open(pwout).read().splitlines()
    if "JOB DONE" not in "\n".join(lines[-5:]) and not any("JOB DONE" in l for l in lines):
        raise SystemExit(f"{pwout}: no JOB DONE -- phonon job incomplete")
    starts = [i for i, l in enumerate(lines) if "Forces acting on atoms" in l]
    if not starts:
        raise SystemExit(f"{pwout}: no forces block")
    forces = []
    for l in lines[starts[-1]:]:
        if "Total force" in l:           # summary line (also contains "force =") -> stop
            break
        if "force =" in l:               # per-atom force line
            p = l.split("force =")[1].split()
            forces.append([float(p[0]), float(p[1]), float(p[2])])
    return np.array(forces) * RY_BOHR_TO_EV_ANG


def parse_volume(pwout):
    """Cell volume (Ang^3) from the QE 'unit-cell volume' line (fixed by cell_dofree='shape')."""
    v = None
    for line in open(pwout):
        if "unit-cell volume" in line:
            v = float(line.split("=")[1].split()[0])   # first token after '=' is the a.u.^3 number
    if v is None:
        raise SystemExit(f"{pwout}: no unit-cell volume")
    return v * BOHR3_TO_ANG3


def free_energy(phase, vidx):
    work = os.path.join(HERE, "runs", phase, f"v{vidx:02d}")
    disp = os.path.join(work, "disp")
    ph = phonopy.load(os.path.join(disp, "phonopy_disp.yaml"))
    n = len(ph.supercells_with_displacements)
    forces = [parse_forces(os.path.join(disp, f"supercell-{i:03d}.out")) for i in range(1, n + 1)]
    ph.forces = forces
    ph.produce_force_constants()
    ph.run_mesh([20, 20, 20])
    ph.run_thermal_properties(t_min=T_MIN, t_max=T_MAX, t_step=T_STEP)
    f_vib = ph.get_thermal_properties_dict()["free_energy"] * KJMOL_TO_EV   # eV/cell, on T-grid
    return parse_volume(os.path.join(work, "static", "pw.out")), f_vib


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", required=True)
    a = ap.parse_args()
    allv = volumes_for(a.phase)
    done = []
    for v in allv:
        try:
            vol, fv = free_energy(a.phase, v)
        except SystemExit as e:
            print(f"[fit_qha] {a.phase} v{v:02d}: SKIP -- {e}")
            continue
        done.append(v)
        print(f"[fit_qha] {a.phase} v{v:02d}: V={vol:.2f} A^3, F_vib(300K)={fv[0]:.4f} eV/cell")
    print(f"[fit_qha] {a.phase}: F_vib built for {len(done)}/{len(allv)} volumes. "
          "EoS fit + eos_params_qha.json emit are TODO.")


if __name__ == "__main__":
    main()
