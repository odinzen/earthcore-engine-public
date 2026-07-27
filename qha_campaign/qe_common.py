"""
Shared Quantum ESPRESSO settings + helpers for the QHA campaign.

Converge ECUTWFC / ECUTRHO / KSPACING / DEGAUSS (and per-element magnetism) in the hcp-Fe
pre-flight, then FREEZE them so every phase is on the same footing. Pseudopotentials are read
from $ESPRESSO_PSEUDO (a dir of .UPF files, e.g. the free SSSP PBE efficiency set); pw.x also
reads that env var, so we don't set pseudo_dir in the input.
"""
import glob
import os
import sys

import numpy as np
from ase import Atoms
from ase.io import write as ase_write

BOHR = 0.529177210903

PSEUDO = os.environ.get("ESPRESSO_PSEUDO")

# --- frozen-after-preflight knobs ---
ECUTWFC = 60.0     # Ry   (TODO converge in pre-flight)
ECUTRHO = 480.0    # Ry   (8x ecutwfc; SSSP US/PAW often want 8-12x)
KSPACING = 0.20    # 1/Ang, 2*pi convention (like VASP KSPACING)
DEGAUSS = 0.02     # Ry (~0.27 eV); Methfessel-Paxton smearing for metals
SMEARING = "mp"
# starting magnetic moment per element (Fe-bearing phases are magnetic at low compression)
MAGMOM = {"Fe": 2.2, "Ni": 0.6, "Si": 0.0, "O": 0.0, "S": 0.0, "C": 0.0}


def pseudos_for(symbols):
    if not PSEUDO:
        sys.exit("ESPRESSO_PSEUDO not set (point it at a dir of .UPF files)")
    out = {}
    for el in sorted(set(symbols)):
        # match files whose name starts with the element symbol then a non-letter (Fe.xxx, Fe_xxx)
        cands = [os.path.basename(f) for f in glob.glob(os.path.join(PSEUDO, "*"))
                 if f.lower().endswith(".upf")]
        hit = [f for f in cands
               if f.upper().startswith(el.upper())
               and (len(f) == len(el) + 4 or not f[len(el)].isalpha())]
        if not hit:
            sys.exit(f"no .UPF for element {el} in {PSEUDO}")
        out[el] = sorted(hit, key=len)[0]
    return out


def kpts_for(atoms, scale=(1, 1, 1)):
    """Monkhorst grid at KSPACING; `scale` divides the grid (use the supercell dim)."""
    recip = 2 * np.pi * np.array(atoms.cell.reciprocal())
    n = [int(np.ceil(np.linalg.norm(recip[i]) / KSPACING)) for i in range(3)]
    return [max(1, n[i] // int(scale[i])) for i in range(3)]


def read_relaxed(pwout):
    """ASE Atoms from the final coordinates of a QE vc-relax output -- bypasses ASE's strict
    read_espresso_out (which asserts on k-points that verbosity=low can omit)."""
    lines = open(pwout).read().splitlines()
    alat = None
    for l in lines:
        if "lattice parameter (alat)" in l:
            alat = float(l.split("=")[1].split()[0]) * BOHR     # bohr -> Ang
    cp = max(i for i, l in enumerate(lines) if l.strip().startswith("CELL_PARAMETERS"))
    u = lines[cp].lower()
    fac = BOHR if "bohr" in u else (alat if "alat" in u and alat else 1.0)
    cell = np.array([[float(x) * fac for x in lines[cp + 1 + k].split()[:3]] for k in range(3)])
    ap = max(i for i, l in enumerate(lines) if l.strip().startswith("ATOMIC_POSITIONS"))
    crystal = "crystal" in lines[ap].lower()
    syms, pos = [], []
    for l in lines[ap + 1:]:
        p = l.split()
        if len(p) < 4 or not p[0][0].isalpha():
            break
        syms.append(p[0]); pos.append([float(p[1]), float(p[2]), float(p[3])])
    return (Atoms(symbols=syms, cell=cell, scaled_positions=pos, pbc=True) if crystal
            else Atoms(symbols=syms, cell=cell, positions=pos, pbc=True))


def write_pw(atoms, dest, calculation, kpts=None, conv_thr=1e-8, nspin=2):
    syms = atoms.get_chemical_symbols()
    a = atoms.copy()
    if nspin == 2:
        a.set_initial_magnetic_moments([MAGMOM.get(s, 0.0) for s in syms])
    data = {
        "control": {"calculation": calculation, "tprnfor": True, "tstress": True,
                    "prefix": "pwscf", "disk_io": "none", "verbosity": "low"},
        "system": {"ecutwfc": ECUTWFC, "ecutrho": ECUTRHO, "occupations": "smearing",
                   "smearing": SMEARING, "degauss": DEGAUSS, "nspin": nspin},
        "electrons": {"conv_thr": conv_thr, "mixing_beta": 0.2, "electron_maxstep": 300,
                      "mixing_mode": "local-TF"},  # stable for magnetic Fe-O/S/C supercells
    }
    if calculation in ("relax", "vc-relax"):
        data["ions"] = {"ion_dynamics": "bfgs"}
    if calculation == "vc-relax":
        data["cell"] = {"cell_dofree": "shape", "press_conv_thr": 0.5}  # fix V, relax shape
    ase_write(dest, a, format="espresso-in", input_data=data,
              pseudopotentials=pseudos_for(syms), kpts=tuple(kpts or kpts_for(a)))
