"""Harvest the DFT method metadata the paper's Methods section needs (no compute, read-only).

Reads the static pw.out headers and the pseudopotential files: QE version, cutoffs, smearing,
per-phase k-grid + relaxed lattice parameter, and the pseudopotential file names + functional/type
(the citable provenance). Run on Sol:  python harvest_meta.py
"""
import glob
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
PSEUDO = os.environ.get("ESPRESSO_PSEUDO", os.path.expanduser("~/pseudo/earthcore"))
PHASES = ["hcp_Fe", "B2_FeSi", "B1_FeO", "Fe3S", "Fe7C3"]


def grep1(path, pat):
    try:
        for line in open(path):
            if re.search(pat, line):
                return line.strip()
    except OSError:
        pass
    return "(not found)"


def main():
    pw = os.path.join(HERE, "runs", "hcp_Fe", "v00", "static", "pw.out")
    print("== code versions ==")
    print("QE:      ", grep1(pw, r"Program PWSCF"))
    try:
        import phonopy
        print("phonopy: ", phonopy.__version__)
    except Exception as e:  # noqa: BLE001
        print("phonopy:  (import failed)", e)

    print("\n== global SCF settings (hcp_Fe v00) ==")
    for label, pat in [("ecutwfc", r"kinetic-energy cutoff"), ("ecutrho", r"charge density cutoff"),
                       ("smearing", r"smearing"), ("degauss", r"gaussian broad")]:
        print(f"  {label:9s}: {grep1(pw, pat)}")

    print("\n== per-phase k-grid + relaxed lattice (static v00) ==")
    for p in PHASES:
        f = os.path.join(HERE, "runs", p, "v00", "static", "pw.out")
        print(f"  {p:8s} {grep1(f, r'number of k points')}  |  {grep1(f, r'lattice parameter')}")

    print(f"\n== pseudopotentials in {PSEUDO} ==")
    files = sorted(glob.glob(os.path.join(PSEUDO, "*.UPF")) + glob.glob(os.path.join(PSEUDO, "*.upf")))
    for f in files:
        head = open(f, errors="ignore").read(6000)
        fx = re.search(r'functional\s*=\s*"?\s*([\w+.-]+)', head, re.I)
        paw = "PAW" if re.search(r'is_paw\s*=\s*"?[Tt]', head) else (
            "USPP" if re.search(r'is_ultrasoft\s*=\s*"?[Tt]', head) else "NC/?")
        print(f"  {os.path.basename(f):40s} functional={fx.group(1) if fx else '?':8s} {paw}")
    if not files:
        print("  (none found -- check ESPRESSO_PSEUDO)")


if __name__ == "__main__":
    main()
