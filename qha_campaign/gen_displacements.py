"""
Stage between the static and phonon arrays (Quantum ESPRESSO). From each relaxed static cell,
generate the phonopy finite-displacement supercells and write each as a QE scf input.

    python gen_displacements.py [--phase hcp_Fe]

Uses phonopy's PYTHON API (not the QE command-line interface, which is fussy about input layout):
build the displaced supercells in memory, write each via ASE as supercell-NNN.in, and save
phonopy_disp.yaml for the force-collection step (fit_qha.py). ASE handles all QE I/O, so units and
format stay consistent (eV/Ang throughout). phonon_jobs.csv is rebuilt from disk each run, so you
can run this per-phase as phases finish.

Run AFTER the static array, BEFORE the phonon array. Submit the phonon array with
--manifest phonon_jobs.csv. Needs ESPRESSO_PSEUDO set (write_pw stamps the pseudo names).
"""
import argparse
import csv
import glob
import os

import numpy as np
from ase import Atoms
from ase.io import read
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms

import qe_common as qe

HERE = os.path.dirname(os.path.abspath(__file__))
# Phases that are nonmagnetic at core compression -> nspin=1 (correct AND stable; nspin=2 supercells
# magnetically diverge for hcp-Fe). The Fe-O/S/C phases stay spin-polarized.
NONMAG = {"hcp_Fe", "B2_FeSi"}


def jobs():
    seen = {}
    with open(os.path.join(HERE, "job_manifest.csv")) as f:
        for r in csv.DictReader(f):
            seen[(r["phase"], int(r["vol_idx"]))] = r["supercell"]   # "4x4x3"
    return seen


def ase_to_phonopy(atoms):
    return PhonopyAtoms(symbols=atoms.get_chemical_symbols(),
                        cell=atoms.cell.array,
                        scaled_positions=atoms.get_scaled_positions())


def phonopy_to_ase(p):
    return Atoms(symbols=p.symbols, cell=p.cell, scaled_positions=p.scaled_positions, pbc=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default=None)
    a = ap.parse_args()
    runroot = os.path.join(HERE, "runs")
    made, skipped = 0, 0
    for (phase, vidx), sc in sorted(jobs().items()):
        if a.phase and phase != a.phase:
            continue
        work = os.path.join(runroot, phase, f"v{vidx:02d}")
        pwout = os.path.join(work, "static", "pw.out")
        if not os.path.isfile(pwout):
            print(f"  skip {phase} v{vidx:02d}: no static pw.out yet"); skipped += 1; continue
        disp = os.path.join(work, "disp"); os.makedirs(disp, exist_ok=True)
        try:
            atoms = read(pwout, format="espresso-out", index=-1)     # relaxed cell, Ang
        except Exception:
            atoms = qe.read_relaxed(pwout)                           # ASE's strict reader choked
        dim = [int(x) for x in sc.split("x")]
        ph = Phonopy(ase_to_phonopy(atoms), supercell_matrix=np.diag(dim))
        ph.generate_displacements(distance=0.01)
        ph.save(os.path.join(disp, "phonopy_disp.yaml"))
        sc_kpts = qe.kpts_for(atoms, scale=dim)        # unit-cell k-density folded for the supercell
        nspin = 1 if phase in NONMAG else 2
        for i, scell in enumerate(ph.supercells_with_displacements, 1):
            qe.write_pw(phonopy_to_ase(scell), os.path.join(disp, f"supercell-{i:03d}.in"),
                        "scf", kpts=sc_kpts, conv_thr=1e-7, nspin=nspin)
        print(f"  {phase} v{vidx:02d}: {len(ph.supercells_with_displacements)} displacements (dim {dim})")
        made += 1

    # rebuild phonon_jobs.csv from disk (incremental-safe)
    rows = []
    for disp in sorted(glob.glob(os.path.join(runroot, "*", "v*", "disp"))):
        parts = disp.split(os.sep); ph_name, vi = parts[-3], int(parts[-2][1:])
        for f in sorted(os.listdir(disp)):
            if f.startswith("supercell-") and f.endswith(".in"):
                rows.append({"phase": ph_name, "vol_idx": vi, "job_kind": "phonon_disp",
                             "disp_idx": int(f.split("-")[1].split(".")[0])})
    if rows:
        with open(os.path.join(HERE, "phonon_jobs.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["phase", "vol_idx", "job_kind", "disp_idx"])
            w.writeheader(); w.writerows(rows)
        print(f"phonon_jobs.csv: {len(rows)} jobs across {len({r['phase'] for r in rows})} phase(s) "
              f"(submit with --manifest phonon_jobs.csv)")
    print(f"done: {made} (phase,vol) sets, {skipped} waiting on static pw.out")


if __name__ == "__main__":
    main()
