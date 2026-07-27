"""Partition double-constraint figures (figs 5,6,7,8). Strict greyscale house style."""
import warnings
warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator

import icb_partition_model as M

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

SOURCE = ("Source: ICB partition model with nickel (icb_partition_model); "
          "PREM Dziewonski & Anderson 1981; partition D_i Alfe 2002, Li 2019, "
          "Badro 2014; shear Murphy 2013, Ono 2013, Prescher 2015. "
          "Greyscale per house figure standard.")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.edgecolor": "black",
    "axes.linewidth": 0.8,
    "axes.grid": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})


def box(ax):
    for s in ("top", "bottom", "left", "right"):
        ax.spines[s].set_visible(True)
        ax.spines[s].set_linewidth(0.8)
        ax.spines[s].set_color("black")
    ax.tick_params(direction="in", color="black", top=True, right=True)


def source_note(fig):
    fig.text(0.5, 0.012, SOURCE, ha="center", va="bottom", fontsize=7,
             color="black", wrap=True)


# Best fit composition
LIQ_O, LIQ_S, LIQ_SI, LIQ_C = 0.18, 0.015, 0.0, 0.0
SOL_O, LIQ_NI = 0.0034, 0.05


# ---------------------------------------------------------------- fig 5
def fig5():
    P = np.linspace(135.0, 329.0, 80)
    prem = np.array([M.prem_outer(p) for p in P])           # rho, vphi
    liq = np.array([M.liquid_props(LIQ_C, LIQ_O, LIQ_SI, LIQ_S, p) for p in P])
    sol = M.solid_props(LIQ_C, LIQ_SI, SOL_O, 0.011, M.P_ICB * 1e9, M.T_ICB)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 4.4))

    # PANEL A density
    axA.plot(P, prem[:, 0], color="black", lw=2.4, ls="-",
             label="PREM outer core")
    axA.plot(P, liq[:, 0], color="0.35", lw=1.6, ls="--",
             label="Best fit liquid, O 18, S 1.5 mol %")
    axA.plot([M.P_ICB], [M.PREM_IC["rho"]], marker="o", ms=9, mfc="white",
             mec="black", mew=1.6, ls="none", label="PREM inner core")
    axA.plot([M.P_ICB], [sol[0]], marker="*", ms=15, mfc="0.25", mec="black",
             mew=1.0, ls="none", label="Partitioned solid")
    # ICB density jump annotation
    jump = 100.0 * (sol[0] - liq[-1, 0]) / liq[-1, 0]
    axA.annotate("", xy=(M.P_ICB, sol[0]), xytext=(M.P_ICB, liq[-1, 0]),
                 arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    axA.text(M.P_ICB - 6, 0.5 * (sol[0] + liq[-1, 0]),
             "ICB density jump\n%.1f %%" % jump, ha="right", va="center",
             fontsize=8, color="black")
    axA.set_xlabel("Pressure, GPa")
    axA.set_ylabel("Density, g cm$^{-3}$")
    axA.set_title("A  DENSITY", loc="left", fontweight="bold", fontsize=9)
    axA.legend(frameon=True, fontsize=7, loc="upper left")
    axA.get_legend().get_frame().set_edgecolor("black")
    axA.get_legend().get_frame().set_linewidth(0.6)
    box(axA)

    # PANEL B bulk sound velocity
    axB.plot(P, prem[:, 1], color="black", lw=2.4, ls="-",
             label="PREM outer core")
    axB.plot(P, liq[:, 1], color="0.35", lw=1.6, ls="--",
             label="Best fit liquid")
    axB.plot([M.P_ICB], [M.PREM_IC["vphi"]], marker="o", ms=9, mfc="white",
             mec="black", mew=1.6, ls="none", label="PREM inner core")
    axB.plot([M.P_ICB], [sol[2]], marker="*", ms=15, mfc="0.25", mec="black",
             mew=1.0, ls="none", label="Partitioned solid")
    axB.set_xlabel("Pressure, GPa")
    axB.set_ylabel("Bulk sound velocity, km s$^{-1}$")
    axB.set_title("B  BULK SOUND VELOCITY", loc="left", fontweight="bold",
                  fontsize=9)
    axB.legend(frameon=True, fontsize=7, loc="upper left")
    axB.get_legend().get_frame().set_edgecolor("black")
    axB.get_legend().get_frame().set_linewidth(0.6)
    box(axB)

    fig.suptitle("Double constraint across the inner core boundary",
                 fontweight="bold", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    source_note(fig)
    out = os.path.join(FIGDIR, "fig5_double_constraint.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


# ---------------------------------------------------------------- fig 6
def fig6():
    # outer-core liquid vs inner-core solid, mole %
    liq = {"O": 18.0, "Si": 0.0, "S": 1.5, "C": 0.0, "Ni": 5.0}
    psol = M.partition_to_solid({"O": 0.18, "Si": 0.0, "S": 0.015, "C": 0.0},
                                )
    sol = {"O": psol["O"] * 100, "Si": psol["Si"] * 100, "S": psol["S"] * 100,
           "C": psol["C"] * 100, "Ni": 5.0}
    els = ["O", "Si", "S", "C", "Ni"]
    x = np.arange(len(els))
    w = 0.38
    liq_v = [max(liq[e], 1e-6) for e in els]
    sol_v = [max(sol[e], 1e-6) for e in els]

    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    b1 = ax.bar(x - w / 2, liq_v, w, color="0.55", edgecolor="black",
                linewidth=0.8, hatch="////", label="Outer core liquid")
    b2 = ax.bar(x + w / 2, sol_v, w, color="0.85", edgecolor="black",
                linewidth=0.8, hatch="....", label="Inner core solid")
    ax.set_yscale("log")
    ax.set_ylim(1e-4, 60)
    ax.set_xticks(x)
    ax.set_xticklabels(els)
    ax.set_xlabel("Light element")
    ax.set_ylabel("Concentration, mole %")
    ax.set_title("A  RESERVOIR CONCENTRATIONS", loc="left", fontweight="bold",
                 fontsize=9)

    # label tiny solid values for O and C
    for e in ("O", "C"):
        i = els.index(e)
        ax.text(i + w / 2, max(sol[e], 1e-6) * 1.4,
                "%.4f" % sol[e] if sol[e] >= 1e-4 else "<1e-4",
                ha="center", va="bottom", fontsize=7, color="black")
    ax.legend(frameon=True, fontsize=8, loc="upper right")
    ax.get_legend().get_frame().set_edgecolor("black")
    ax.get_legend().get_frame().set_linewidth(0.6)
    box(ax)
    fig.suptitle("Partition of light elements, liquid versus solid",
                 fontweight="bold", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    source_note(fig)
    out = os.path.join(FIGDIR, "fig6_partition_bars.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


# ---------------------------------------------------------------- fig 7
def fig7():
    pc = M.PP["partition_coefficients"]
    els = ["C", "O", "Si", "S", "Ni"]
    D = [pc[e]["D"] for e in els]
    shades = ["0.2", "0.4", "0.55", "0.7", "0.85"]
    hatches = ["////", "\\\\\\\\", "xxxx", "....", "++++"]
    y = np.arange(len(els))

    fig, ax = plt.subplots(figsize=(8.0, 4.4))
    for i, (yy, d) in enumerate(zip(y, D)):
        ax.barh(yy, d, color=shades[i], edgecolor="black", linewidth=0.8,
                hatch=hatches[i])
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(els)
    ax.set_xlabel("Partition coefficient D$_i$ = x$_{solid}$ / x$_{liquid}$")
    ax.set_ylabel("Light element")
    ax.axvline(1.0, color="black", lw=1.4, ls="-")
    ax.set_xlim(1e-5, 1e1)
    # labels left/right of D=1
    ax.text(2.0, len(els) - 0.4, "favors solid", ha="left", va="center",
            fontsize=8, color="black", style="italic")
    ax.text(0.5, len(els) - 0.4, "favors liquid", ha="right", va="center",
            fontsize=8, color="black", style="italic")
    for i, d in enumerate(D):
        ax.text(d * 1.5, i, "%.2g" % d, va="center", ha="left", fontsize=7,
                color="black")
    ax.set_title("A  PARTITION COEFFICIENTS", loc="left", fontweight="bold",
                 fontsize=9)
    box(ax)
    fig.suptitle("Freezing front partition coefficients",
                 fontweight="bold", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    source_note(fig)
    out = os.path.join(FIGDIR, "fig7_partition_coeffs.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


# ---------------------------------------------------------------- fig 8
def fig8():
    df = pd.read_csv(os.path.join(HERE, "icb_partition_results.csv"))
    sl = df[(df.liq_x_si == 0) & (df.liq_x_c == 0)].copy()
    sl["o_pct"] = sl.liq_x_o * 100.0
    sl["s_pct"] = sl.liq_x_s * 100.0

    ovals = np.sort(sl.o_pct.unique())
    svals = np.sort(sl.s_pct.unique())
    Z = np.full((len(svals), len(ovals)), np.nan)
    oidx = {v: i for i, v in enumerate(ovals)}
    sidx = {v: i for i, v in enumerate(svals)}
    for _, r in sl.iterrows():
        Z[sidx[r.s_pct], oidx[r.o_pct]] = r.rms_misfit_pct
    Og, Sg = np.meshgrid(ovals, svals)

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    cf = ax.contourf(Og, Sg, Z, levels=14, cmap="Greys")
    adm = sl[sl.admissible]
    ax.plot(adm.o_pct, adm.s_pct, marker="o", ls="none", mfc="black",
            mec="black", ms=6, label="Admissible")
    ax.plot([18.0], [1.5], marker="*", ms=18, mfc="white", mec="black",
            mew=1.6, ls="none", label="Best fit")
    ax.set_xlabel("Liquid oxygen, mole %")
    ax.set_ylabel("Liquid sulfur, mole %")
    ax.set_title("A  ADMISSIBLE REGION", loc="left", fontweight="bold",
                 fontsize=9)
    leg = ax.legend(frameon=True, fontsize=8, loc="upper right")
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.6)
    box(ax)
    cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("RMS misfit, %", fontsize=8)
    cb.outline.set_linewidth(0.6)
    cb.outline.set_edgecolor("black")
    cb.ax.tick_params(width=0.6, color="black")

    fig.suptitle("PREM misfit over liquid oxygen and sulfur",
                 fontweight="bold", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    source_note(fig)
    out = os.path.join(FIGDIR, "fig8_admissible_region.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


if __name__ == "__main__":
    for fn in (fig5, fig6, fig7, fig8):
        print("wrote", fn())
