"""Diagnose why B1_FeO's QHA thermal is broken (negative Grueneisen).

Reprocesses the EXISTING B1_FeO force sets (no new DFT): reports the minimum phonon frequency per
volume (a negative value = imaginary/soft mode -> QHA invalid, a structural instability, not a
convergence bug) and compares F_vib(300 K) from the original shifted [20,20,20] mesh against a
gamma-centered mesh (to see if the mesh shift alone is the culprit).

    python diag_feo.py
"""
import os

import numpy as np
import phonopy

import fit_qha as F

PHASE = "B1_FeO"


def main():
    print(f"{'vol':4s} {'V(A^3)':>7s} {'min_freq(THz)':>13s} {'Fvib300_shift':>14s} "
          f"{'Fvib300_gamma':>14s}  note")
    for v in F.volumes_for(PHASE):
        work = os.path.join(F.HERE, "runs", PHASE, f"v{v:02d}")
        disp = os.path.join(work, "disp")
        ph = phonopy.load(os.path.join(disp, "phonopy_disp.yaml"))
        n = len(ph.supercells_with_displacements)
        ph.forces = [F.parse_forces(os.path.join(disp, f"supercell-{i:03d}.out"))
                     for i in range(1, n + 1)]
        ph.produce_force_constants()

        ph.run_mesh([24, 24, 24], is_gamma_center=True)
        fmin = float(np.min(ph.get_mesh_dict()["frequencies"]))
        ph.run_thermal_properties(t_min=300, t_max=300, t_step=100)
        fg = ph.get_thermal_properties_dict()["free_energy"][0] * F.KJMOL_TO_EV

        ph.run_mesh([20, 20, 20])
        ph.run_thermal_properties(t_min=300, t_max=300, t_step=100)
        fs = ph.get_thermal_properties_dict()["free_energy"][0] * F.KJMOL_TO_EV

        vol = F.parse_volume(os.path.join(work, "static", "pw.out"))
        note = "IMAGINARY MODE" if fmin < -0.1 else ("soft" if fmin < 0.1 else "")
        print(f"v{v:02d}  {vol:7.2f} {fmin:13.3f} {fs:14.4f} {fg:14.4f}  {note}")
    print("\nnegative min_freq at several volumes -> B1_FeO is dynamically unstable there, so its "
          "QHA F(V,T) is not physical (use literature thermal or B8_FeO). If frequencies are all "
          "positive and F_vib_gamma is smoother than shift, the mesh was the problem -> reprocess.")


if __name__ == "__main__":
    main()
