"""
Validation gate at CORE conditions: compare QE's pressure at each sampled (compressed) volume to
the phase's sourced equation of state in ../core/eos_params.json (Vinet from V0,K0,K'). Run after a
phase's static array finishes.

    python fit_ev.py --phase hcp_Fe     # Dewaele 2006
    python fit_ev.py --phase B2_FeSi    # Fischer 2014
    python fit_ev.py --phase B8_FeO     # Fischer 2011

Our volume grid is V/V0 = 0.55-0.78 (core pressures), so the zero-pressure V0/K0 are OUTSIDE the
data -- fitting them would extrapolate. Instead we check P(V) directly: QE vs the literature Vinet
at each sampled volume. Agreement within a few percent across the core range validates the
pseudopotential / cutoff / k-grid for that phase.
"""
import argparse
import glob
import json
import math
import os

from ase.io import read

HERE = os.path.dirname(os.path.abspath(__file__))
NA = 6.02214076e23
A3_TO_CM3MOL = 1e-24 * NA            # Ang^3/atom -> cm3/mol-atom


def phase_eos(phase):
    em = json.load(open(os.path.join(HERE, "..", "core", "eos_params.json")))["endmembers"]
    if phase not in em:
        raise SystemExit(f"{phase} not in eos_params.json (Fe3S lives in partition_params.json)")
    e = em[phase]
    return (e["V_0"]["value"], e["K_0"]["value"], e["Kprime_0"]["value"],
            e["molar_mass"]["n_atoms"], e["V_0"].get("source", "?"))


def vinet_P(vv0, K0, Kp):
    x = vv0 ** (1.0 / 3.0)
    return 3.0 * K0 * (1.0 - x) / x ** 2 * math.exp(1.5 * (Kp - 1.0) * (1.0 - x))


def qe_pressure_gpa(pwout):
    p = None
    with open(pwout) as f:
        for line in f:
            if "P=" in line:                       # QE total-stress line, kbar
                try:
                    p = float(line.split("P=")[1].split()[0])
                except (IndexError, ValueError):
                    pass
    return None if p is None else p / 10.0          # kbar -> GPa


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", default="hcp_Fe")
    a = ap.parse_args()
    V0, K0, Kp, n_form, src = phase_eos(a.phase)
    V0_atom = V0 / n_form                            # cm3/mol-atom
    outs = sorted(glob.glob(os.path.join(HERE, "runs", a.phase, "v*", "static", "pw.out")))
    if not outs:
        raise SystemExit(f"no static pw.out for {a.phase}")
    print(f"\n=== {a.phase}: QE pressure vs {src} (Vinet K0={K0} K'={Kp}) over {len(outs)} vols ===")
    print(f"  {'V/V0':>6} {'P_QE':>9} {'P_lit':>9} {'diff':>8}")
    diffs = []
    for o in outs:
        atoms = read(o, format="espresso-out", index=-1)
        v_atom = (atoms.get_volume() / len(atoms)) * A3_TO_CM3MOL
        pqe = qe_pressure_gpa(o)
        if pqe is None:
            print(f"  {v_atom/V0_atom:6.3f}   (no pressure -- did it finish?)"); continue
        plit = vinet_P(v_atom / V0_atom, K0, Kp)
        d = 100.0 * (pqe - plit) / plit
        diffs.append(abs(d))
        print(f"  {v_atom/V0_atom:6.3f} {pqe:9.1f} {plit:9.1f} {d:+7.1f}%")
    if diffs:
        mad = sum(diffs) / len(diffs)
        print(f"\n  {len(diffs)}/{len(outs)} finished; mean |diff| = {mad:.1f}%  ->  "
              + ("PASS (validated at core conditions)" if mad < 8
                 else "CHECK settings / pseudopotential"))


if __name__ == "__main__":
    main()
