"""
Build the QHA campaign job matrix for the FeSiOCS core endmembers.

Reads V_0/K_0 per phase from ../core/eos_params.json (the paper's
sourced ledger) and combines them with per-phase STRUCTURAL facts (cell size,
supercell, displacement count) to emit the full (phase, volume, displacement) job
list plus a Sol-budget summary. Pure stdlib; no DFT, no numpy. Run anywhere.

    python make_manifest.py

Writes job_manifest.csv and job_manifest.json next to this file and prints the budget.

NOTE: the resource numbers are PLANNING ESTIMATES. Replace them with the measured
pre-flight values (hcp-Fe, one volume) before trusting the campaign-wide total --
see PLAN.md "Pre-flight".
"""
import csv
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
LEDGER = os.path.join(HERE, "..", "core", "eos_params.json")

# Volume grid: compression ratios V/V0. v00-v08 (0.55-0.78) bracket ~100-400 GPa (CMB 136 ->
# center 364) and carry the phonons for F(V,T). v09-v11 (0.86-1.02) are the low-compression
# ANCHOR added 2026-07-11: STATIC E(V) only (no phonons -- the HP phases FeSi/FeO can go
# dynamically unstable near ambient, but their static energy is well-defined), extending the cold
# curve to P~0 so ambient V0/K0/K0' are pinned instead of extrapolated from compressed-only data.
V_RATIOS = [0.55, 0.58, 0.61, 0.64, 0.67, 0.70, 0.73, 0.76, 0.78, 0.86, 0.94, 1.02]
ANCHOR_FROM_IDX = 9        # v09+ are static-only anchor volumes; do NOT gen phonons for these

# T-grid is POST-PROCESSING (phonon DOS -> F(V,T)); it costs no extra DFT jobs.
T_GRID = {"min_K": 300, "max_K": 7000, "step_K": 100,
          "note": "covers the 5500-6200 K ICB scan plus margin; free (no DFT per T)."}

# Per-phase structural facts. atoms_supercell / n_disp drive cost. CONFIRM these against
# the actual relaxed structures before the production run -- they set the memory ceiling.
# n_disp = symmetry-inequivalent finite displacements phonopy will generate (estimate).
PHASES = {
    "hcp_Fe": {"ledger_key": "hcp_Fe", "space_group": "P6_3/mmc", "atoms_prim": 2,
               "supercell": "4x4x3", "atoms_supercell": 96, "n_disp": 1,
               "role": "iron base / BENCHMARK vs Dewaele 2006", "tier": "small"},
    "B2_FeSi": {"ledger_key": "B2_FeSi", "space_group": "Pm-3m", "atoms_prim": 2,
                "supercell": "4x4x4", "atoms_supercell": 128, "n_disp": 1,
                "role": "silicon arm", "tier": "small"},
    "B8_FeO": {"ledger_key": "B8_FeO", "space_group": "P6_3/mmc (NiAs)", "atoms_prim": 4,
               "supercell": "3x3x2", "atoms_supercell": 72, "n_disp": 3,
               "role": "OXYGEN ARM (load-bearing) -- default FeO", "tier": "small"},
    "B1_FeO": {"ledger_key": "B1_FeO", "space_group": "Fm-3m (NaCl)", "atoms_prim": 2,
               "supercell": "3x3x3", "atoms_supercell": 54, "n_disp": 1,
               "role": "oxygen arm (stable-phase alt, Fischer 2011)", "tier": "small"},
    "Fe3S": {"ledger_key": None, "space_group": "Pnma (MP mp-1189032)", "atoms_prim": 16,
             "supercell": "2x2x2", "atoms_supercell": 128, "n_disp": 24,
             "role": "sulfur arm (params live in partition_params.json)", "tier": "medium"},
    "Fe7C3": {"ledger_key": "Fe7C3", "space_group": "P6_3mc (MP mp-18257)", "atoms_prim": 20,
              "supercell": "2x2x2", "atoms_supercell": 160, "n_disp": 12,
              "role": "carbon arm -- HIGHEST-UNCERTAINTY phase (eos FLAG); biggest win, GATE IT",
              "tier": "large"},
}

# Crude resource model vs supercell atom count. Tune from the pre-flight.
# cost scales ~ N^3 for DFT; anchor at a small cell and a per-job wall budget.
def estimate(atoms_supercell):
    n = atoms_supercell
    cores = 32 if n <= 96 else 48            # one node; metals like dense k-grids
    core_hours = 0.02 * (n / 32.0) ** 3 * 32 # rough; anchored ~0.6 core-h at 32 atoms
    wall_h = core_hours / cores
    mem_gb = 0.9 + 0.05 * n                  # rough peak per job
    return cores, round(core_hours, 2), round(wall_h, 2), round(mem_gb, 1)


def load_v0(key):
    if key is None:
        return None
    with open(LEDGER) as f:
        em = json.load(f)["endmembers"]
    return em[key]["V_0"]["value"]  # cm3/mol (atom or formula; ratio grid is unitless)


def build():
    jobs = []
    for phase, sp in PHASES.items():
        v0 = load_v0(sp["ledger_key"])
        cores, ch, wall, mem = estimate(sp["atoms_supercell"])
        for vi, r in enumerate(V_RATIOS):
            # job 0 = the static E-V SCF; jobs 1..n_disp = phonon displacements.
            # anchor volumes (v09+) are static-only -- no phonon displacements.
            n_disp = 0 if vi >= ANCHOR_FROM_IDX else sp["n_disp"]
            for d in range(n_disp + 1):
                kind = "static_EV" if d == 0 else "phonon_disp"
                jobs.append({
                    "phase": phase, "tier": sp["tier"], "role": sp["role"],
                    "supercell": sp["supercell"], "atoms_supercell": sp["atoms_supercell"],
                    "vol_idx": vi, "V_over_V0": r,
                    "V0_cm3_mol": v0, "V_cm3_mol": (round(v0 * r, 4) if v0 else None),
                    "job_kind": kind, "disp_idx": (None if d == 0 else d),
                    "cores": cores, "est_core_hours": ch, "est_wall_h": wall, "est_mem_gb": mem,
                    "over_wall_budget": wall > 8.0,
                })
    return jobs


def summarize(jobs):
    by_phase = {}
    for j in jobs:
        p = by_phase.setdefault(j["phase"], {"jobs": 0, "core_hours": 0.0, "max_wall": 0.0,
                                             "max_mem": 0.0, "over": 0})
        p["jobs"] += 1
        p["core_hours"] += j["est_core_hours"]
        p["max_wall"] = max(p["max_wall"], j["est_wall_h"])
        p["max_mem"] = max(p["max_mem"], j["est_mem_gb"])
        p["over"] += int(j["over_wall_budget"])
    return by_phase


if __name__ == "__main__":
    jobs = build()
    with open(os.path.join(HERE, "job_manifest.json"), "w") as f:
        json.dump({"v_ratios": V_RATIOS, "t_grid": T_GRID, "phases": PHASES, "jobs": jobs}, f, indent=2)
    with open(os.path.join(HERE, "job_manifest.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(jobs[0].keys()))
        w.writeheader()
        w.writerows(jobs)

    by_phase = summarize(jobs)
    print(f"Volume grid: {len(V_RATIOS)} points, V/V0 in [{V_RATIOS[0]}, {V_RATIOS[-1]}]")
    print(f"T-grid (free, post-proc): {T_GRID['min_K']}-{T_GRID['max_K']} K\n")
    print(f"{'phase':10s} {'tier':7s} {'jobs':>5s} {'core-h':>8s} {'maxWall(h)':>11s} {'maxMem(GB)':>11s} {'>8h':>4s}")
    tot_jobs = tot_ch = 0
    for p, s in by_phase.items():
        print(f"{p:10s} {PHASES[p]['tier']:7s} {s['jobs']:5d} {s['core_hours']:8.1f} "
              f"{s['max_wall']:11.2f} {s['max_mem']:11.1f} {s['over']:4d}")
        tot_jobs += s["jobs"]
        tot_ch += s["core_hours"]
    print(f"{'TOTAL':10s} {'':7s} {tot_jobs:5d} {tot_ch:8.1f}")
    over = [j for j in jobs if j["over_wall_budget"]]
    if over:
        phs = sorted({j['phase'] for j in over})
        print(f"\nWARNING: {len(over)} jobs exceed the 8 h wall budget (phases: {', '.join(phs)}).")
        print("  -> shrink that supercell or split the displacement array before scaling out (PLAN.md Phase A).")
    print("\nGATE: Fe7C3 is the highest-uncertainty phase AND the largest cell -- run its pre-flight last,")
    print("      and only after hcp-Fe validates against Dewaele 2006.")
    print("Wrote job_manifest.csv and job_manifest.json. Numbers are ESTIMATES until the pre-flight.")
