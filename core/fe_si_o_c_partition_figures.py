"""Partition double-constraint figures (figs 5-10) from production/ (full_space_rerun.py --production).
Strict greyscale house style."""
import warnings
warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
import json

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


# Best fit composition: one source for every figure, the velocities and the text.
PROD = os.path.join(HERE, "production")
BEST = json.load(open(os.path.join(PROD, "best_fit.json")))
LIQ_O, LIQ_S, LIQ_SI, LIQ_C = (BEST["liquid_x"][k] for k in ("O", "S", "Si", "C"))
SOL_O, SOL_S, SOL_SI, SOL_C = (BEST["solid_x"][k] for k in ("O", "S", "Si", "C"))
MW = {"Fe": 55.845, "Ni": 58.693, "O": 15.999, "Si": 28.0855, "S": 32.06, "C": 12.011}


def wt_pct(xo, xsi, xs, xc, x_ni=0.05):
    x = {"O": xo, "Si": xsi, "S": xs, "C": xc, "Ni": x_ni}
    x["Fe"] = 1.0 - sum(x.values())
    tot = sum(x[k] * MW[k] for k in x)
    return {k: 100 * x[k] * MW[k] / tot for k in x}


def legend(ax, **kw):
    leg = ax.legend(frameon=True, fontsize=8, **kw)
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.6)


# ---------------------------------------------------------------- fig 5
def fig5():
    P = np.linspace(135.0, 329.0, 80)
    # The real PREM profile (make_prem_outer_profile.py); the fit itself uses only its CMB and ICB values.
    prem_tab = np.loadtxt(os.path.join(HERE, "prem_outer_core_profile.csv"), delimiter=",")
    P_prem, prem = prem_tab[:, 0], prem_tab[:, 1:]          # rho, vphi
    liq = np.array([M.liquid_props(LIQ_C, LIQ_O, LIQ_SI, LIQ_S, p) for p in P])
    sol = M.solid_props(SOL_C, SOL_SI, SOL_O, SOL_S, M.P_ICB * 1e9, M.T_ICB)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 4.4))

    # PANEL A density
    axA.plot(P_prem, prem[:, 0], color="black", lw=2.4, ls="-",
             label="PREM outer core")
    axA.plot(P, liq[:, 0], color="0.35", lw=1.6, ls="--",
             label="Best fit liquid")
    axA.plot([M.P_ICB], [M.PREM_IC["rho"]], marker="o", ms=9, mfc="white",
             mec="black", mew=1.6, ls="none", label="PREM inner core")
    axA.plot([M.P_ICB], [sol[0]], marker="*", ms=15, mfc="0.25", mec="black",
             mew=1.0, ls="none", label="Partitioned solid")
    axA.set_xlabel("Pressure, GPa")
    axA.set_ylabel("Density, g cm$^{-3}$")
    axA.set_title("(a)", loc="left", fontweight="bold", fontsize=9)
    axA.legend(frameon=True, fontsize=7, loc="upper left")
    axA.get_legend().get_frame().set_edgecolor("black")
    axA.get_legend().get_frame().set_linewidth(0.6)
    box(axA)

    # PANEL B bulk sound velocity
    axB.plot(P_prem, prem[:, 1], color="black", lw=2.4, ls="-",
             label="PREM outer core")
    axB.plot(P, liq[:, 1], color="0.35", lw=1.6, ls="--",
             label="Best fit liquid")
    axB.plot([M.P_ICB], [M.PREM_IC["vphi"]], marker="o", ms=9, mfc="white",
             mec="black", mew=1.6, ls="none", label="PREM inner core")
    axB.plot([M.P_ICB], [sol[2]], marker="*", ms=15, mfc="0.25", mec="black",
             mew=1.0, ls="none", label="Partitioned solid")
    axB.set_xlabel("Pressure, GPa")
    axB.set_ylabel("Bulk sound velocity, km s$^{-1}$")
    axB.set_title("(b)", loc="left", fontweight="bold",
                  fontsize=9)
    axB.legend(frameon=True, fontsize=7, loc="upper left")
    axB.get_legend().get_frame().set_edgecolor("black")
    axB.get_legend().get_frame().set_linewidth(0.6)
    box(axB)

    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    out = os.path.join(FIGDIR, "fig5_double_constraint.png")
    fig.savefig(out, dpi=600, metadata={"Software": None})
    plt.close(fig)
    return out


# ---------------------------------------------------------------- fig 6
def fig6():
    # outer-core liquid vs inner-core solid, wt%
    liq = {**BEST["liquid_wt_pct"], "Ni": wt_pct(LIQ_O, LIQ_SI, LIQ_S, LIQ_C)["Ni"]}
    sol = {**BEST["solid_wt_pct"], "Ni": wt_pct(SOL_O, SOL_SI, SOL_S, SOL_C)["Ni"]}
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
    ax.set_ylim(1e-2, 20)
    ax.set_xticks(x)
    ax.set_xticklabels(els)
    ax.set_xlabel("Element")
    ax.set_ylabel("Concentration, wt%")

    ax.legend(frameon=True, fontsize=8, loc="upper right")
    ax.get_legend().get_frame().set_edgecolor("black")
    ax.get_legend().get_frame().set_linewidth(0.6)
    box(ax)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    out = os.path.join(FIGDIR, "fig6_partition_bars.png")
    fig.savefig(out, dpi=600, metadata={"Software": None})
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
    ax.set_ylabel("Element")
    ax.axvline(1.0, color="black", lw=1.4, ls="-")
    ax.set_xlim(1e-5, 1e1)
    # region labels either side of D = 1, kept inside the frame
    ax.set_ylim(-0.6, len(els) + 0.1)
    ax.text(2.0, len(els) - 0.3, "favors solid", ha="left", va="center",
            fontsize=8, color="black")
    ax.text(0.5, len(els) - 0.3, "favors liquid", ha="right", va="center",
            fontsize=8, color="black")
    box(ax)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    out = os.path.join(FIGDIR, "fig7_partition_coeffs.png")
    fig.savefig(out, dpi=600, metadata={"Software": None})
    plt.close(fig)
    return out


# ---------------------------------------------------------------- fig 8
def fig8():
    """Oxygen against carbon at the best-fit Si and S with the carbon cap lifted: where the cap sits."""
    import full_space_rerun as F
    F._M = M
    F._CAPS_WT = {}
    M.base._phase_props.cache_clear()
    grid = dict(O=np.round(np.arange(0, 0.205, 0.01), 3), Si=np.array([LIQ_SI]),
                S=np.array([LIQ_S]), C=np.round(np.arange(0, 0.205, 0.01), 3))
    df = F._sweep(grid, M.T_ICB)
    w = [wt_pct(r.liq_x_o, r.liq_x_si, r.liq_x_s, r.liq_x_c) for r in df.itertuples()]
    df["wO"], df["wC"] = [x["O"] for x in w], [x["C"] for x in w]
    Zo = df.pivot(index="liq_x_c", columns="liq_x_o", values="wO").values
    Zc = df.pivot(index="liq_x_c", columns="liq_x_o", values="wC").values
    Zr = df.pivot(index="liq_x_c", columns="liq_x_o", values="rms_misfit_pct").values

    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    cf = ax.contourf(Zo, Zc, Zr, levels=14, cmap="Greys")
    a = df[df.admissible]
    ax.plot(a.wO, a.wC, marker="o", ls="none", mfc="black", mec="black", ms=5, label="Admissible")
    cap = BEST["caps_wt_pct"]["C"]
    ax.axhline(cap, color="black", lw=1.2, ls="--", label=f"Carbon limit, {cap:g} wt%")
    bw = BEST["liquid_wt_pct"]
    ax.plot([bw["O"]], [bw["C"]], marker="*", ms=18, mfc="white", mec="black", mew=1.6, ls="none",
            label="Best fit")
    ax.set_xlim(0, 5)
    ax.set_ylim(0, 3)
    ax.set_xlabel("Liquid oxygen, wt%")
    ax.set_ylabel("Liquid carbon, wt%")
    legend(ax, loc="upper right")
    box(ax)
    cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("RMS misfit, %", fontsize=8)
    cb.outline.set_linewidth(0.6)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "fig8_admissible_region.png")
    fig.savefig(out, dpi=600, metadata={"Software": None})
    plt.close(fig)
    return out


# ---------------------------------------------------------------- fig 9
def fig9():
    """Monte Carlo: best-fit outer core oxygen and the boundary density jump over admissible draws."""
    u = pd.read_csv(os.path.join(PROD, "uq_draws.csv"))
    a = u[u.admissible == True]
    o_wt = np.array([wt_pct(r.best_o, r.best_si, r.best_s, r.best_c)["O"] for r in a.itertuples()])
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 4.0))
    axA.hist(o_wt, bins=np.arange(0, 5.51, 0.25), color="0.65", edgecolor="black", lw=0.6)
    for q, ls in ((50, "-"), (16, "--"), (84, "--")):
        axA.axvline(np.percentile(o_wt, q), color="black", lw=1.2, ls=ls)
    axA.set_xlabel("Best fit outer core oxygen, wt%")
    axA.set_ylabel("Admissible draws")
    axA.set_title("(a)", loc="left", fontweight="bold", fontsize=9)
    box(axA)
    axB.hist(a.best_jump, bins=np.arange(4.0, 8.51, 0.25), color="0.65", edgecolor="black", lw=0.6)
    axB.axvspan(4.5, 4.8, color="0.85", zorder=0)
    axB.axvline(np.percentile(a.best_jump, 50), color="black", lw=1.2)
    axB.set_xlabel("Boundary density jump, %")
    axB.set_ylabel("Admissible draws")
    axB.set_title("(b)", loc="left", fontweight="bold", fontsize=9)
    box(axB)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "fig9_uq.png")
    fig.savefig(out, dpi=600, metadata={"Software": None})
    plt.close(fig)
    return out


# ---------------------------------------------------------------- fig 10
SCEN = [("T_ICB_5500K", "5500 K"), ("T_ICB_5700K", "5700 K"), ("baseline_5800K", "5800 K"),
        ("T_ICB_6000K", "6000 K"), ("T_ICB_6200K", "6200 K"), ("B1_oxide", "B1 oxide"),
        ("D_Si_0.67", "D$_{Si}$ 0.67"), ("D_Si_1.00", "D$_{Si}$ 1.00")]


def fig10():
    """Sensitivity: admissible outer core oxygen range and best fit per scenario."""
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    for i, (key, lab) in enumerate(SCEN):
        d = pd.read_csv(os.path.join(PROD, f"sweep_{key}.csv"))
        a = d[d.admissible]
        o = [wt_pct(r.liq_x_o, r.liq_x_si, r.liq_x_s, r.liq_x_c)["O"] for r in a.itertuples()]
        b = a.loc[a.rms_misfit_pct.idxmin()]
        ax.plot([min(o), max(o)], [i, i], color="0.45", lw=6, solid_capstyle="butt")
        ax.plot([wt_pct(b.liq_x_o, b.liq_x_si, b.liq_x_s, b.liq_x_c)["O"]], [i], marker="*", ms=13,
                mfc="white", mec="black", mew=1.2, ls="none")
    ax.set_yticks(range(len(SCEN)))
    ax.set_yticklabels([lab for _, lab in SCEN])
    ax.invert_yaxis()
    ax.set_xlim(0, 5)
    ax.set_xlabel("Outer core oxygen, wt%")
    ax.plot([], [], color="0.45", lw=6, label="Admissible range")
    ax.plot([], [], marker="*", ms=13, mfc="white", mec="black", ls="none", label="Best fit")
    legend(ax, loc="upper left")
    box(ax)
    fig.tight_layout()
    out = os.path.join(FIGDIR, "fig10_sensitivity.png")
    fig.savefig(out, dpi=600, metadata={"Software": None})
    plt.close(fig)
    return out


if __name__ == "__main__":
    for fn in (fig5, fig6, fig7, fig8, fig9, fig10):
        print("wrote", fn())
