"""Fig 11: aggregate Vp, Vphi, Vs and Poisson, best fit solid vs PREM. Greyscale."""
import warnings
warnings.filterwarnings("ignore")
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import icb_velocity as V

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


def main():
    r = V.vp_vs(0.0, 0.0, 0.0034, 0.011)   # x_c, x_si, x_o, x_s
    prem = V.SC["PREM_inner_core"]
    prem_vphi = 10.26

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(9.2, 4.4),
                                   gridspec_kw={"width_ratios": [3, 1.3]})

    # PANEL A: Vp, Vphi, Vs
    labels = ["Vp", "V$\\phi$", "Vs"]
    model = [r["vp"], r["vphi"], r["vs"]]
    premv = [prem["Vp_km_s"], prem_vphi, prem["Vs_km_s"]]
    x = np.arange(len(labels))
    w = 0.38
    axA.bar(x - w / 2, model, w, color="0.55", edgecolor="black",
            linewidth=0.8, hatch="////", label="Best fit solid")
    axA.bar(x + w / 2, premv, w, color="0.85", edgecolor="black",
            linewidth=0.8, hatch="....", label="PREM inner core")
    axA.set_xticks(x)
    axA.set_xticklabels(labels)
    axA.set_ylabel("Velocity, km s$^{-1}$")
    axA.set_xlabel("Wave mode")
    axA.set_title("A  AGGREGATE VELOCITIES", loc="left", fontweight="bold",
                  fontsize=9)
    # annotate the Vs overshoot
    dvs = 100.0 * (r["vs"] - prem["Vs_km_s"]) / prem["Vs_km_s"]
    axA.text(2, max(model[2], premv[2]) + 0.4, "Vs %+.0f %%" % dvs,
             ha="center", va="bottom", fontsize=8, color="black")
    axA.set_ylim(0, 13)
    leg = axA.legend(frameon=True, fontsize=8, loc="upper right")
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.6)
    box(axA)

    # PANEL B: Poisson
    xb = np.arange(2)
    pv = [r["poisson"], prem["poisson"]]
    axB.bar([0], [pv[0]], 0.6, color="0.55", edgecolor="black", linewidth=0.8,
            hatch="////", label="Best fit solid")
    axB.bar([1], [pv[1]], 0.6, color="0.85", edgecolor="black", linewidth=0.8,
            hatch="....", label="PREM inner core")
    axB.set_xticks(xb)
    axB.set_xticklabels(["Model", "PREM"])
    axB.set_ylabel("Poisson ratio")
    axB.set_xlabel("Reservoir")
    axB.set_ylim(0, 0.5)
    axB.set_title("B  POISSON RATIO", loc="left", fontweight="bold",
                  fontsize=9)
    for i, v in enumerate(pv):
        axB.text(i, v + 0.01, "%.2f" % v, ha="center", va="bottom",
                 fontsize=8, color="black")
    box(axB)

    fig.suptitle("Shear anomaly of the inner core solid",
                 fontweight="bold", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.text(0.5, 0.012, SOURCE, ha="center", va="bottom", fontsize=7,
             color="black", wrap=True)
    out = os.path.join(FIGDIR, "fig11_vp_vs_poisson.png")
    fig.savefig(out, dpi=200)
    plt.close(fig)
    print("wrote", out)
    print("model vp=%.2f vphi=%.2f vs=%.2f poisson=%.2f" %
          (r["vp"], r["vphi"], r["vs"], r["poisson"]))


if __name__ == "__main__":
    main()
