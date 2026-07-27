"""
Per-job driver (Quantum ESPRESSO). One SLURM array index -> one pw.x run.

    python run_one.py --index $SLURM_ARRAY_TASK_ID --manifest job_manifest.csv

  static_EV   : scale the starter cell to V/V0, vc-relax at FIXED volume (cell_dofree='shape').
  phonon_disp : scf+forces on a phonopy-generated supercell input (gen_displacements.py).

Runs ON Sol inside the qha conda env (ase + phonopy + qe). Env:
  QE_CMD           launch command (default 'srun pw.x')
  ESPRESSO_PSEUDO  dir of .UPF pseudopotentials (e.g. SSSP PBE efficiency)
Outputs under runs/<phase>/v<idx>/. Idempotent via DONE markers.
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys

from ase.io import read

import qe_common as qe

HERE = os.path.dirname(os.path.abspath(__file__))
QE_CMD = os.environ.get("QE_CMD", "srun pw.x").split()
POSCAR = {p: os.path.join(HERE, "structures", f"POSCAR_{p}")
          for p in ("hcp_Fe", "B2_FeSi", "B1_FeO", "B8_FeO", "Fe3S", "Fe7C3")}


def row_at(manifest, index):
    with open(manifest) as f:
        rows = list(csv.DictReader(f))
    if not 0 <= index < len(rows):
        sys.exit(f"index {index} out of range (manifest has {len(rows)} jobs)")
    return rows[index]


def run_pw(workdir, infile, outfile, marker):
    if os.path.exists(os.path.join(workdir, marker)):
        print(f"[run_one] already done: {workdir}/{marker}"); return
    print(f"[run_one] pw.x in {workdir}: {' '.join(QE_CMD)} -in {infile}", flush=True)
    with open(os.path.join(workdir, outfile), "w") as out:
        r = subprocess.run(QE_CMD + ["-in", infile], cwd=workdir, stdout=out)
    ok = r.returncode == 0 and _job_done(os.path.join(workdir, outfile))
    if not ok:
        sys.exit(f"pw.x failed (rc={r.returncode}) in {workdir} -- see {outfile}")
    open(os.path.join(workdir, marker), "w").close()


def _job_done(outpath):
    try:
        with open(outpath) as f:
            txt = f.read()
        # require BOTH job completion AND a forces block -- a non-converged SCF prints JOB DONE
        # but no forces, and those outputs are useless for phonons/relaxation.
        return "JOB DONE" in txt and "Forces acting on atoms" in txt
    except OSError:
        return False


def do_static(row, work):
    base = POSCAR[row["phase"]]
    if not os.path.isfile(base):
        sys.exit(f"no starter POSCAR for {row['phase']} ({base}); see structures/README.md")
    wd = os.path.join(work, "static"); os.makedirs(wd, exist_ok=True)
    atoms = read(base, format="vasp")
    f = float(row["V_over_V0"]) ** (1.0 / 3.0)
    atoms.set_cell(atoms.cell * f, scale_atoms=True)
    qe.write_pw(atoms, os.path.join(wd, "pw.in"), "vc-relax")
    run_pw(wd, "pw.in", "pw.out", "DONE")
    print(f"[run_one] static done; relaxed cell + E(V) in {wd}/pw.out")


def do_phonon(row, work):
    d = int(row["disp_idx"])
    disp = os.path.join(work, "disp")
    inf = f"supercell-{d:03d}.in"
    if not os.path.isfile(os.path.join(disp, inf)):
        sys.exit(f"{disp}/{inf} missing -- run gen_displacements.py for "
                 f"{row['phase']} v{row['vol_idx']} after its static pw.out exists")
    run_pw(disp, inf, f"supercell-{d:03d}.out", f".done-{d:03d}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--index", type=int, required=True)
    ap.add_argument("--manifest", default=os.path.join(HERE, "job_manifest.csv"))
    ap.add_argument("--runroot", default=os.path.join(HERE, "runs"))
    a = ap.parse_args()
    row = row_at(a.manifest, a.index)
    work = os.path.join(a.runroot, row["phase"], f"v{int(row['vol_idx']):02d}")
    os.makedirs(work, exist_ok=True)
    print(f"[run_one] index={a.index} phase={row['phase']} v{row['vol_idx']} "
          f"V/V0={row.get('V_over_V0', '-')} kind={row['job_kind']} disp={row['disp_idx']}", flush=True)
    (do_static if row["job_kind"] == "static_EV" else do_phonon)(row, work)


if __name__ == "__main__":
    main()
