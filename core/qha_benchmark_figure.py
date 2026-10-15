"""Figure: first-principles validation of the solid equation of state (FeSiOCS paper, Fig 12).

Left panel  -- calculated hcp-Fe density vs pressure at 300 K, against the Dewaele 2006 EoS and PREM.
Right panel -- density difference (calculated minus imported) vs pressure for each solid phase.

Greyscale, 600 DPI, no in-figure title, labels defined in the caption (Kristina figure standard).
Reads the QHA tables in ../qha_campaign/. Run: python qha_benchmark_figure.py
"""
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
QHA = os.path.join(HERE, "..", "qha_campaign")
EVA3_GPA = 160.2176621
NA = 6.02214076e23
C = 1e24 / NA
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

# phase: (mass_cell g/mol, lit V0 cell A^3, lit K0 GPa, lit K0', line style, grey)
UPC = {"hcp_Fe": 2, "B2_FeSi": 1, "B1_FeO": 1, "Fe3S": 4, "Fe7C3": 2}
LIT = {
    "hcp_Fe": (2 * 55.845, 6.753 * C * 2, 163.4, 5.38, "-", "0.0", "hcp Fe"),
    "B2_FeSi": (55.845 + 28.0855, 12.828 * C, 230.6, 4.17, "--", "0.3", "B2 FeSi"),
    "B1_FeO": (55.845 + 15.999, 12.256 * C, 149.4, 3.60, "-.", "0.45", "B1 FeO"),
    "Fe3S": (4 * (3 * 55.845 + 32.06), 56.76 * C * 2, 156.0, 3.80, ":", "0.55", "Fe3S"),
    "Fe7C3": (2 * (7 * 55.845 + 3 * 12.011), 53.0 * C * 2, 291.0, 4.50, (0, (3, 1, 1, 1)), "0.35", "Fe7C3"),
}

plt.rcParams.update({"font.size": 13, "axes.labelsize": 14, "legend.fontsize": 10,
                     "xtick.labelsize": 12, "ytick.labelsize": 12})


def vinet_P(V, V0, B0, B1):
    x = (V / V0) ** (1 / 3)
    return 3 * B0 * (1 - x) / x**2 * np.exp(1.5 * (B1 - 1) * (1 - x))


def _box(ax):
    for s in ax.spines.values():
        s.set_color("black")
        s.set_linewidth(0.8)
    ax.tick_params(direction="out", length=4, width=0.7, colors="black")
    ax.grid(False)


def load(ph):
    d = json.load(open(os.path.join(QHA, f"qha_{ph}.json")))
    vs = sorted(d["volumes"], key=lambda z: z["V_ang3"])
    V = np.array([z["V_ang3"] for z in vs])
    iT = d["temperatures_K"].index(300.0)
    F = np.array([z["E_static_eV"] + z["F_vib_eV"][iT] for z in vs])
    return V, F


fig, (axL, axR) = plt.subplots(1, 2, figsize=(11.4, 5.0))

# left: hcp density vs P, calculated vs Dewaele vs PREM
mass, lV0, lK0, lKp, _, _, _ = LIT["hcp_Fe"]
V, F = load("hcp_Fe")
P = np.array([-(F[i + 1] - F[i - 1]) / (V[i + 1] - V[i - 1]) * EVA3_GPA for i in range(1, len(V) - 1)])
rho_qha = mass / NA / (V[1:-1] * 1e-24)
Pg = np.linspace(80, 360, 120)
rho_dew = np.array([mass / NA / (brentq(lambda v: vinet_P(v, lV0, lK0 / EVA3_GPA, lKp) * EVA3_GPA - p,
                                        lV0 * 0.4, lV0 * 1.05) * 1e-24) for p in Pg])
prem = 13.0885 - 8.8381 * (363.85 - Pg) / 952.0
axL.plot(Pg, rho_dew, ls="--", color="0.45", lw=2.2, label="Dewaele et al. (2006)")
axL.plot(P, rho_qha, ls="none", marker="o", ms=7, mfc="white", mec="black", mew=1.4,
         label="calculated (present work)")
axL.plot(Pg, prem, ls=":", color="0.0", lw=1.4, label="PREM inner core")
axL.set_xlabel("Pressure (GPa)")
axL.set_ylabel("Density of hcp iron (g/cm$^3$)")
axL.set_xlim(80, 360)
_box(axL)
axL.legend(loc="upper left", frameon=True, edgecolor="black")

# right: density difference vs P per phase
for ph, (mass, lV0, lK0, lKp, ls, grey, lab) in LIT.items():
    V, F = load(ph)
    diffs, Ps = [], []
    for i in range(1, len(V) - 1):
        p = -(F[i + 1] - F[i - 1]) / (V[i + 1] - V[i - 1]) * EVA3_GPA
        if p <= 1:
            continue
        Vl = brentq(lambda v: vinet_P(v, lV0, lK0 / EVA3_GPA, lKp) * EVA3_GPA - p, lV0 * 0.35, lV0 * 1.2)
        diffs.append(100 * (Vl - V[i]) / V[i])  # (rho_qha-rho_lit)/rho_lit = (V_lit-V_qha)/V_qha
        Ps.append(p)
    axR.plot(Ps, diffs, ls=ls, color=grey, lw=2.2, label=lab)
axR.axhline(0.0, color="0.0", lw=0.8)
axR.set_xlabel("Pressure (GPa)")
axR.set_ylabel("Density difference, calculated $-$ imported (%)")
axR.set_xlim(60, 480)
_box(axR)
axR.legend(loc="lower right", frameon=True, edgecolor="black", ncol=1)

fig.tight_layout()
fig.savefig(os.path.join(FIGDIR, "fig12_qha_benchmark.png"), dpi=600, bbox_inches="tight", metadata={"Software": None})
fig.savefig(os.path.join(FIGDIR, "fig12_qha_benchmark.pdf"), bbox_inches="tight", metadata={"Creator": None, "Producer": None})
print("wrote fig12_qha_benchmark.{png,pdf}")
