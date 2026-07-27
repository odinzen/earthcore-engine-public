"""Greyscale figure suite for the oxygen-provenance (magma-ocean) paper.

Kristina figure standard: no in-figure titles, no annotation boxes or arrows, no
source blocks. The caption carries every number. Black and white only, 600 DPI,
axis labels legible without zooming. Series are distinguished by line style, grey
shade, hatching, and legend entries. All numbers come from oxygen_solver.

Run:  python magma_figures.py
"""
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import oxygen_solver as S

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

PB = 54.0
TB = 3350.0
FEO_B = S.P["mantle_FeO"]["present_wt_pct"]   # present-day BSE, 8.05 wt%
DPI = 600

plt.rcParams.update({
    "font.size": 13,
    "axes.labelsize": 14,
    "legend.fontsize": 11,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
})

SHADE = {"a": "0.15", "b": "0.45", "c": "0.62"}
LS = {"a": "-", "b": "--", "c": "-."}


def _box(ax):
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color("black")
        s.set_linewidth(0.8)
    ax.tick_params(direction="out", length=4, width=0.7, colors="black")
    ax.set_facecolor("white")
    ax.grid(False)


def _legend(ax, **kw):
    leg = ax.legend(frameon=True, **kw)
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.8)
    leg.get_frame().set_facecolor("white")
    return leg


def _save(fig, stem):
    fig.savefig(os.path.join(FIGDIR, stem + ".png"), dpi=DPI, bbox_inches="tight")
    fig.savefig(os.path.join(FIGDIR, stem + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def fig1_two_routes():
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ox = S.solve_oxidized_route(5.0)
    feo_ox = ox["feo"]
    Ts = np.linspace(3000.0, 4600.0, 60)
    for feo, key, lab in ((FEO_B, "a", "mantle FeO 8 wt%"),
                          (feo_ox, "b", f"mantle FeO {feo_ox:.0f} wt%")):
        wt = [S.core_O_wt(T, PB, feo, mode="table2") for T in Ts]
        ax.plot(Ts, wt, color=SHADE[key], ls=LS[key], lw=2.4, label=lab)

    ax.axhline(5.0, color="0.0", lw=1.0, ls=(0, (1, 1)), label="core target, 5 wt% O")
    o_ss = S.core_O_wt(TB, PB, FEO_B, mode="table2")
    ax.plot(TB, o_ss, marker="o", ms=8, mfc="white", mec="black", mew=1.4, ls="none",
            label="single-stage equilibration")
    t_hot = S.solve_T_for_O(5.0, PB, FEO_B, mode="table2")
    ax.plot(t_hot, 5.0, marker="s", ms=8, mfc="0.15", mec="black", mew=1.0, ls="none",
            label="5 wt% O, hot route (forbidden)")
    ax.plot(ox["T"], 5.0, marker="^", ms=9, mfc="0.45", mec="black", mew=1.0, ls="none",
            label="5 wt% O, oxidized route (on liquidus)")

    ax.set_xlabel("Equilibration temperature (K)")
    ax.set_ylabel("Oxygen in core-forming metal (wt%)")
    ax.set_xlim(3000, 4600)
    ax.set_ylim(0, 7)
    _box(ax)
    _legend(ax, loc="lower right")
    _save(fig, "fig1_two_routes")


def fig2_fo2():
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    feos = np.linspace(5.0, 30.0, 40)
    ideal, frost = [], []
    for feo in feos:
        xfe = S.solve_metal(TB, PB, feo, mode="table2")["X"]["Fe"]
        ideal.append(S.delta_IW(feo, xfe))
        frost.append(S.delta_IW(feo, xfe, T=TB, Pg=PB, activity=True))
    ax.plot(feos, ideal, color=SHADE["b"], ls=LS["b"], lw=2.4, label="mole-fraction approximation")
    ax.plot(feos, frost, color=SHADE["a"], ls=LS["a"], lw=2.4, label="non-ideal activities")
    ax.axhline(-2.0, color="0.0", lw=1.0, ls=(0, (1, 1)), label="iron-wüstite minus 2")

    ax.set_xlabel("Mantle FeO (wt%)")
    ax.set_ylabel("Oxygen fugacity, delta-IW (log units)")
    ax.set_xlim(5, 30)
    _box(ax)
    _legend(ax, loc="lower right")
    _save(fig, "fig2_fo2")


def fig3_uq():
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    feos = np.linspace(6.0, 28.0, 12)
    med, lo, hi = [], [], []
    for feo in feos:
        u = S.uq_O_at_state(TB, PB, feo, n=1500)
        med.append(u["median"])
        lo.append(u["lo"])
        hi.append(u["hi"])
    ax.fill_between(feos, lo, hi, facecolor="0.82", edgecolor="0.4", lw=0.6,
                    hatch="////", label="16th to 84th percentile")
    ax.plot(feos, med, color="0.0", ls="-", lw=2.2, label="median")
    ax.errorbar(FEO_B, 1.6, yerr=0.3, fmt="o", ms=8, mfc="white", mec="black",
                mew=1.4, ecolor="black", capsize=5, ls="none",
                label="single-stage equilibration")

    ax.set_xlabel("Mantle FeO (wt%)")
    ax.set_ylabel("Oxygen in core-forming metal (wt%)")
    ax.set_xlim(6, 28)
    ax.set_ylim(0, None)
    _box(ax)
    _legend(ax, loc="upper left")
    _save(fig, "fig3_uq")


def fig4_mantle_magma():
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    o = np.linspace(1.0, 8.0, 40)
    init = np.array([S.mantle_FeO_before_core(x)[0] for x in o])
    ax.plot(o, init, color="0.0", ls="-", lw=2.4, label="pre-core mantle FeO")
    ax.axhspan(20.0, 25.0, facecolor="0.85", edgecolor="0.4", lw=0.6, hatch="\\\\",
               label="oxidized early mantle")
    ax.axhline(S.P["mantle_FeO"]["present_wt_pct"], color="0.45", ls="-.", lw=1.6,
               label="present-day mantle FeO")
    i5 = S.mantle_FeO_before_core(5.0)[0]
    ax.plot(5.0, i5, marker="o", ms=8, mfc="white", mec="black", mew=1.4, ls="none",
            label="5 wt% O core")

    ax.set_xlabel("Oxygen in the core (wt%)")
    ax.set_ylabel("Pre-core mantle FeO (wt%)")
    ax.set_xlim(1, 8)
    ax.set_ylim(6, 30)
    _box(ax)
    _legend(ax, loc="lower right")
    _save(fig, "fig4_mantle_magma")


def fig5_reconciliation():
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    f = np.linspace(0.45, 1.0, 160)
    feo = np.array([S.early_mantle_FeO(5.0, fi) for fi in f])
    ax.plot(f, feo, color="0.0", ls="-", lw=2.4, label="5 wt% O core, mass balance")

    ax.axhline(28.0, color="0.45", ls="--", lw=1.8, label="oxidized-route requirement")

    ax.plot(1.0, S.early_mantle_FeO(5.0, 1.0), marker="s", ms=8, mfc="white", mec="black",
            mew=1.4, ls="none", label="whole-mantle value")
    fstar = S.acquisition_fraction(5.0, 28.0)
    ax.plot(fstar, 28.0, marker="o", ms=9, mfc="0.15", mec="black", mew=1.0, ls="none",
            label="acquisition fraction")

    ax.set_xlabel("Mantle mass at oxygen acquisition (fraction of present)")
    ax.set_ylabel("Local silicate FeO (wt%)")
    ax.set_xlim(0.45, 1.0)
    ax.set_ylim(15, 45)
    _box(ax)
    _legend(ax, loc="upper right")
    _save(fig, "fig5_reconciliation")


if __name__ == "__main__":
    fig1_two_routes()
    fig2_fo2()
    fig3_uq()
    fig4_mantle_magma()
    fig5_reconciliation()
    print("wrote 5 figures (png + pdf) to", FIGDIR)
